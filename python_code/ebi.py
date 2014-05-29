#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT
from os.path import basename, exists, join
from os import mkdir
from datetime import date, timedelta
from xml.etree import ElementTree as ET
from xml.dom import minidom
from xml.sax.saxutils import escape

from make_per_sample_fastq import split_helper


class InvalidMetadataError(Exception):
    pass


class SampleAlreadyExistsError(Exception):
    pass


class NoXMLError(Exception):
    pass


def clean_whitespace(s):
    return ' '.join(s.split())


def iter_file_via_list_of_dicts(input_file):
    """Iterates over a TSV file, yielding dicts keyed by the column headers

    Parameters
    ----------
    input_file : file
        The file to iterate over

    Returns
    -------
    generator
        Yields dicts keyed by the column headers

    Notes
    -----
    - Assumes the first line contains the column headers
    - Assumes no duplicate column headers
    - All column headers will be transformed to lowercase
    - Skips lines containing only whitespace
    """
    headers = [x.strip().lower() for x in input_file.readline().split('\t')]
    for line in input_file:
        line = line.strip()
        if not line:
            continue

        line_elements = line.split('\t')
        yield dict(zip(headers, line_elements))


class EBISubmission(object):
    def __init__(self, study_id, study_title, study_abstract,
                 investigation_type, empty_value='unknown', **kwargs):
        self.study_id = study_id
        self.study_title = study_title
        self.study_abstract = study_abstract
        self.investigation_type = investigation_type

        self.study_xml_fp = None
        self.sample_xml_fp = None
        self.experiment_xml_fp = None
        self.run_xml_fp = None

        # dicts that map investigation_type to library attributes
        lib_strategies = {'metagenome': 'POOLCLONE',
                          'mimarks-survey': 'AMPLICON'}
        lib_selections = {'mimarks-survey': 'PCR'}
        lib_sources = {}

        self.library_strategy = lib_strategies.get(
            self.investigation_type, "OTHER")
        self.library_source = lib_sources.get(
            self.investigation_type, "METAGENOMIC")
        self.library_selection = lib_selections.get(
            self.investigation_type, "unspecififed")

        # This is allows addition of other arbitrary study metadata
        try:
            self.additional_metadata = {
                str(k): str(v) if v is not None else empty_value
                for k, v in kwargs.iteritems()}
        except ValueError:
            raise InvalidMetadataError("All additional metadata passed via "
                                       "kwargs to the EBISubmission "
                                       "constructor must be strings")

        # This will hold the submission's samples, keyed by the sample name
        self.samples = {}

    def _get_study_alias(self):
        """Format alias using ``self.study_id``"""
        return 'qiime_study_' + escape(clean_whitespace(str(self.study_id)))

    def _get_sample_alias(self, sample_name):
        """Format alias using ``self.study_id``, `sample_name`"""
        return "%s:%s" % (self._get_study_alias(),
                          escape(clean_whitespace(str(sample_name))))

    def _get_experiment_alias(self, sample_name, row_number):
        """Format alias using ``self.study_id``, `sample_name`, `row_number`

        `row_number` comes from the index of the prep in the sample's prep
        list.
        """
        return "%s:%d" % (self._get_sample_alias(sample_name),
                          row_number)

    def _get_submission_alias(self):
        """Format alias using ``self.study_id``"""
        safe_study_id = escape(clean_whitespace(str(self.study_id)))
        return 'qiime_submission_' + safe_study_id

    def _get_library_name(self, sample_name, row_number):
        """Format alias using `sample_name`, `row_number`

        `row_number` comes from the index of the prep in the sample's prep
        list.
        """
        return '%s:%d' % (escape(clean_whitespace(sample_name)), row_number)

    def _add_dict_as_tags_and_values(self, parent_node, attribute_element_name,
                                     data_dict):
        """Format key/value data using a common EBI XML motif"""
        for attr, val in sorted(data_dict.iteritems()):
            attribute_element = ET.SubElement(parent_node,
                                              attribute_element_name)
            tag = ET.SubElement(attribute_element, 'TAG')
            tag.text = escape(clean_whitespace(attr))
            value = ET.SubElement(attribute_element, 'VALUE')
            value.text = escape(clean_whitespace(val))

    def generate_study_xml(self):
        """Generates the study XML file

        Returns
        -------
        xml.etree.Element
            The root elelement of the generated ``ElementTree``
        """
        study_set = ET.Element('STUDY_SET', {
            'xmlns:xsi': "http://www.w3.org/2001/XMLSchema-instance",
            'xsi:noNamespaceSchemaLocation': "ftp://ftp.sra.ebi.ac.uk/meta/xsd"
                                             "/sra_1_3/SRA.study.xsd"})

        study = ET.SubElement(study_set, 'STUDY', {
            'alias': self._get_study_alias(),
            'center_name': "CCME-COLORADO"}
        )

        descriptor = ET.SubElement(study, 'DESCRIPTOR')
        study_title = ET.SubElement(descriptor, 'STUDY_TITLE')
        study_title.text = escape(clean_whitespace(self.study_title))
        study_type = ET.SubElement(descriptor, 'STUDY_TYPE', {
            'existing_study_type': escape(clean_whitespace(
                self.investigation_type))}
        )
        study_abstract = ET.SubElement(descriptor, 'STUDY_ABSTRACT')
        study_abstract.text = clean_whitespace(escape(self.study_abstract))

        if self.additional_metadata:
            study_attributes = ET.SubElement(study, 'STUDY_ATTRIBUTES')
            self._add_dict_as_tags_and_values(study_attributes,
                                              'STUDY_ATTRIBUTE',
                                              self.additional_metadata)

        return study_set

    def add_sample(self, sample_name, taxon_id=None, description=None,
                   empty_value='unknown', **kwargs):
        """Adds sample information to the current submission

        Parameters
        ----------
        sample_name : str
            Unique identifier for the sample
        taxon_id : str, optional
            Defaults to ``None``. If not provided, the `empty_value` will be
            used for the taxon ID
        description : str, optional
            Defaults to ``None``. If not provided, the `empty_value` will be
            used for the description
        empty_value : str, optional
            Defaults to "unknown". This value will be used for `taxon_id` or
            `description` when they are not supplied, and for all ``None``
            values in ``kwargs``

        Raises
        ------
        InvalidMetadataError
            If metadata is passed via `kwargs` that is not interpretable as
            a string
        SampleAlreadyExistsError
            If `sample_name` already exists in the ``samples`` dict
        """
        if sample_name in self.samples:
            raise SampleAlreadyExistsError("Two samples with the same "
                                           "sample_name cannot be added to "
                                           "the same submission. "
                                           "(sample_name: %s)" % sample_name)

        self.samples[sample_name] = {}

        self.samples[sample_name]['taxon_id'] = empty_value if \
            taxon_id is None else taxon_id
        self.samples[sample_name]['taxon_id'] = \
            escape(clean_whitespace(self.samples[sample_name]['taxon_id']))

        self.samples[sample_name]['description'] = empty_value if \
            description is None else description
        self.samples[sample_name]['description'] = \
            escape(clean_whitespace(self.samples[sample_name]['description']))

        try:
            self.samples[sample_name]['attributes'] = {
                str(k): str(v) if v is not None else empty_value
                for k, v in kwargs.iteritems()}

        except ValueError:
            raise InvalidMetadataError("All metadata passed to add_sample "
                                       "via kwargs must be strings.")

        self.samples[sample_name]['preps'] = []

    def generate_sample_xml(self):
        """Generates the sample XML file

        Returns
        -------
        xml.etree.Element
            The root elelement of the generated ``ElementTree``
        """
        sample_set = ET.Element('SAMPLE_SET', {
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xsi:noNamespaceSchemaLocation": "ftp://ftp.sra.ebi.ac.uk/meta/xsd"
                                             "/sra_1_3/SRA.sample.xsd"})

        for sample_name, sample_info in sorted(self.samples.iteritems()):
            sample = ET.SubElement(sample_set, 'SAMPLE', {
                'alias': self._get_sample_alias(sample_name),
                'center_name': 'CCME-COLORADO'}
            )

            sample_title = ET.SubElement(sample, 'TITLE')
            sample_title.text = escape(clean_whitespace(sample_name))

            sample_name_element = ET.SubElement(sample, 'SAMPLE_NAME')
            taxon_id = ET.SubElement(sample_name_element, 'TAXON_ID')
            taxon_id.text = escape(clean_whitespace(sample_info['taxon_id']))

            description = ET.SubElement(sample, 'DESCRIPTION')
            description.text = escape(clean_whitespace(
                sample_info['description']))

            if sample_info['attributes']:
                sample_attributes = ET.SubElement(sample, 'SAMPLE_ATTRIBUTES')
                self._add_dict_as_tags_and_values(sample_attributes,
                                                  'SAMPLE_ATTRIBUTE',
                                                  sample_info['attributes'])

        return sample_set

    def add_sample_prep(self, sample_name, platform, file_type, file_path,
                        experiment_design_description,
                        library_construction_protocol,
                        empty_value='unknown', **kwargs):
        """Add prep info for an existing sample

        Parameters
        ----------
        sample_name : str
            The `sample_name` must exist in the ``samples`` list
        platform : {'LS454', 'ILLUMINA', 'UNKNOWN'}
            The sequencing platform
        file_type : {'sff', 'fastq', 'fasta'}
            The type of file
        file_path : str
            The path to the sequence file for this sample
        experiment_design_description : str
        library_construction_protocol : str
        empty_value : str, optional
            Defaults to "unknown". This value will be used for all ``None``
            values in ``kwargs``

        Raises
        ------
        InvalidMetadataError
            If metadata is passed via `kwargs` that is not interpretable as
            a string
        KeyError
            If `sample_name` is not in the list of samples in the
            ``EBISubmission`` object
        """
        try:
            prep_info = {str(k): str(v) if v is not None else empty_value
                         for k, v in kwargs.iteritems()}
        except ValueError:
            raise InvalidMetadataError("All metadata passed to "
                                       "add_sample_prep via kwargs must be "
                                       "strings.")

        prep_info['platform'] = platform
        prep_info['file_type'] = file_type
        prep_info['file_path'] = file_path
        prep_info['experiment_design_description'] = \
            experiment_design_description
        prep_info['library_construction_protocol'] = \
            library_construction_protocol

        self.samples[sample_name]['preps'].append(prep_info)

    def _generate_library_descriptor(self, design, sample_name, row_number,
                                     library_construction_protocol):
        """This XML element (and its subelements) must be written for every
        sample, but its generation depends on only study-level information.
        Therefore, we can breka it out into its own method.
        """

        library_descriptor = ET.SubElement(design, 'LIBRARY_DESCRIPTOR')
        library_name = ET.SubElement(library_descriptor, 'LIBRARY_NAME')
        library_name.text = self._get_library_name(sample_name,
                                                   row_number)
        library_strategy = ET.SubElement(library_descriptor,
                                         "LIBRARY_STRATEGY")
        library_strategy.text = self.library_strategy
        library_source = ET.SubElement(library_descriptor,
                                       "LIBRARY_SOURCE")
        library_source.text = self.library_source
        library_selection = ET.SubElement(library_descriptor,
                                          "LIBRARY_SELECTION")
        library_selection.text = self.library_selection
        library_layout = ET.SubElement(library_descriptor,
                                       "LIBRARY_LAYOUT")
        single = ET.SubElement(library_layout, "SINGLE")
        library_construction_protocol_element = ET.SubElement(
            library_descriptor, "LIBRARY_CONSTRUCTION_PROTOCOL")
        library_construction_protocol_element.text = escape(clean_whitespace(
            library_construction_protocol))

    def _generate_spot_descriptor(self, design, platform):
        """This XML element (and its subelements) must be written for every
        sample, but its generation depends on only study-level information.
        Therefore, we can breka it out into its own method.
        """
        # This section applies only to the LS454 platform
        if platform is not 'LS454':
            return

        # There is some hard-coded information in here, but this is what we
        # have always done in the past...
        spot_descriptor = ET.SubElement(design, 'SPOT_DESCRIPTOR')
        spot_decode_spec = ET.SubElemenet(spot_descriptor, 'SPOT_DECODE_SPEC')
        read_spec = ET.SubElemenet(spot_descriptor, 'READ_SPEC')

        read_index = ET.SubElemenet(read_spec, 'READ_INDEX')
        read_index.text = '0'
        read_class = ET.SubElement(read_spec, 'READ_CLASS')
        read_class.text = 'Application Read'
        read_type = ET.SubElement(read_spec, 'READ_TYPE')
        read_type.text = 'Forward'
        base_coord = ET.SubElement(read_spec, 'BASE_COORD')
        base_coord.text = '1'

    def generate_experiment_xml(self):
        """Generates the experiment XML file

        Returns
        -------
        xml.etree.Element
            The root elelement of the generated ``ElementTree``
        """
        study_alias = self._get_study_alias()
        experiment_set = ET.Element('EXPERIMENT_SET', {
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xsi:noNamespaceSchemaLocation": "ftp://ftp.sra.ebi.ac.uk/meta/xsd"
                                             "/sra_1_3/SRA.experiment.xsd"})
        for sample_name, sample_info in sorted(self.samples.iteritems()):
            sample_alias = self._get_sample_alias(sample_name)
            for row_number, prep_info in enumerate(sample_info['preps']):
                experiment_alias = self._get_experiment_alias(sample_name,
                                                              row_number)
                platform = prep_info['platform']
                experiment = ET.SubElement(experiment_set, 'EXPERIMENT', {
                    'alias': experiment_alias,
                    'center_name': 'CCME-COLORADO'}
                )
                title = ET.SubElement(experiment, 'TITLE')
                title.text = experiment_alias
                study_ref = ET.SubElement(experiment, 'STUDY_REF', {
                    'refname': study_alias}
                )

                design = ET.SubElement(experiment, 'DESIGN')
                design_description = ET.SubElement(design,
                                                   'DESIGN_DESCRIPTION')
                design_description.text = escape(clean_whitespace(
                    prep_info['experiment_design_description']))
                sample_descriptor = ET.SubElement(
                    design, 'SAMPLE_DESCRIPTOR', {'refname': sample_alias}
                )

                self._generate_library_descriptor(
                    design, sample_name, row_number,
                    prep_info['library_construction_protocol']
                )

                self._generate_spot_descriptor(design, platform)

                platform_element = ET.SubElement(experiment, 'PLATFORM')
                platform_info = ET.SubElement(platform_element,
                                              platform.upper())
                instrument_model = ET.SubElement(platform_info,
                                                 'INSTRUMENT_MODEL')
                instrument_model.text = 'unspecified'

                if prep_info:
                    experiment_attributes = ET.SubElement(
                        experiment, 'EXPERIMENT_ATTRIBUTES')
                    self._add_dict_as_tags_and_values(experiment_attributes,
                                                      'EXPERIMENT_ATTRIBUTE',
                                                      prep_info)

        return experiment_set

    def generate_run_xml(self):
        """Generates the run XML file

        Returns
        -------
        xml.etree.Element
            The root elelement of the generated ``ElementTree``
        """
        # TODO: md5 stuff
        run_set = ET.Element('RUN_SET', {
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xsi:noNamespaceSchemaLocation": "ftp://ftp.sra.ebi.ac.uk/meta/xsd"
                                             "/sra_1_3/SRA.run.xsd"})
        for sample_name, sample_info in sorted(self.samples.iteritems()):
            sample_alias = self._get_sample_alias(sample_name)
            for row_number, prep_info in enumerate(sample_info['preps']):
                experiment_alias = self._get_experiment_alias(sample_name,
                                                              row_number)
                platform = prep_info['platform']
                file_type = prep_info['file_type']
                file_path = prep_info['file_path']

                run = ET.SubElement(run_set, 'RUN', {
                    'alias': basename(file_path) + '_run',
                    'center_name': 'CCME-COLORADO'}
                )
                experiment_ref = ET.SubElement(run, 'EXPERIMENT_REF', {
                    'refname': experiment_alias}
                )
                data_block = ET.SubElement(run, 'DATA_BLOCK')
                files = ET.SubElement(data_block, 'FILES')
                file_element = ET.SubElement(files, 'FILE', {
                    'filename': basename(file_path),
                    'filetype': file_type,
                    'quality_scring_system': 'phred',
                    'checksum_method': 'MD5',
                    'checksum': 'NONE'}  # TODO: checksum
                )

        return run_set

    def generate_submission_xml(self, action):
        """Generates the submission XML file

        Parameters
        ----------
        action : {'ADD', 'VALIDATE', 'UPDATE'}
            What action to take when communicating with EBI

        Returns
        -------
        xml.etree.Element
            The root elelement of the generated ``ElementTree``

        Raises
        ------
        NoXMLError
            If one of the necessary XML files has not been generated
        """
        if any([self.study_xml_fp is None,
                self.sample_xml_fp is None,
                self.experiment_xml_fp is None,
                self.run_xml_fp is None]):
            raise NoXMLError("One of the necessary XML files has not been "
                             "generated. Make sure you write out the other "
                             "XML files before attempting to write the "
                             "submission XML file.")

        submission_set = ET.Element('SUBMISSION_SET', {
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xsi:noNamespaceSchemaLocation": "ftp://ftp.sra.ebi.ac.uk/meta/xsd"
                                             "/sra_1_3/SRA.submission.xsd"})
        submission = ET.SubElement(submission_set, 'SUBMISSION', {
            'alias': self._get_submission_alias(),
            'center_name': 'CCME-COLORADO'}
        )
        actions = ET.SubElement(submission, 'ACTIONS')
        study_action = ET.SubElement(actions, 'ACTION', {
            'schema': 'study'}
        )
        sample_action = ET.SubElement(actions, 'ACTION', {
            'schema': 'sample'}
        )
        experiment_action = ET.SubElement(actions, 'ACTION', {
            'schema': 'experiment'}
        )
        run_action = ET.SubElement(actions, 'ACTION', {
            'schema': 'run'}
        )

        if action is 'ADD':
            hold_action = ET.SubElement(actions, 'HOLD', {
                'HoldUntilDate': str(date.today() + timedelta(365))}
            )

        return submission_set

    def write_study_xml(self, fp):
        """Write the study XML file using the current data

        Parameters
        ----------
        fp : str
            The filepath to which the XML will be written

        Notes
        -----
        If `fp` points to an existing file, it will be overwritten
        """
        study_xml_element = self.generate_study_xml()

        study_xml = minidom.parseString(ET.tostring(
            study_xml_element))

        with open(fp, 'w') as outfile:
            outfile.write(study_xml.toprettyxml(indent='  ', encoding='UTF-8'))

        self.study_xml_fp = fp

    def write_sample_xml(self, fp):
        """Write the sample XML file using the current data

        Parameters
        ----------
        fp : str
            The filepath to which the XML will be written

        Notes
        -----
        If `fp` points to an existing file, it will be overwritten
        """
        sample_xml_element = self.generate_sample_xml()

        sample_xml = minidom.parseString(ET.tostring(
            sample_xml_element))

        with open(fp, 'w') as outfile:
            outfile.write(sample_xml.toprettyxml(indent='  ',
                                                 encoding='UTF-8'))

        self.sample_xml_fp = fp

    def write_experiment_xml(self, fp):
        """Write the experiment XML file using the current data

        Parameters
        ----------
        fp : str
            The filepath to which the XML will be written

        Notes
        -----
        If `fp` points to an existing file, it will be overwritten
        """
        experiment_xml_element = self.generate_experiment_xml()

        experiment_xml = minidom.parseString(ET.tostring(
            experiment_xml_element))

        with open(fp, 'w') as outfile:
            outfile.write(experiment_xml.toprettyxml(indent='  ',
                                                     encoding='UTF-8'))

        self.experiment_xml_fp = fp

    def write_run_xml(self, fp):
        """Write the run XML file using the current data

        Parameters
        ----------
        fp : str
            The filepath to which the XML will be written

        Notes
        -----
        If `fp` points to an existing file, it will be overwritten
        """
        run_xml_element = self.generate_run_xml()

        run_xml = minidom.parseString(ET.tostring(
            run_xml_element))

        with open(fp, 'w') as outfile:
            outfile.write(run_xml.toprettyxml(indent='  ', encoding='UTF-8'))

        self.run_xml_fp = fp

    def write_submission_xml(self, fp, action):
        """Write the submission XML file using the current data

        Parameters
        ----------
        fp : str
            The filepath to which the XML will be written
        action : {'ADD', 'VALIDATE', 'UPDATE'}
            What action to take when communicating with EBI

        Notes
        -----
        If `fp` points to an existing file, it will be overwritten
        """
        submission_xml_element = self.generate_submission_xml(action)

        submission_xml = minidom.parseString(ET.tostring(
            submission_xml_element))

        with open(fp, 'w') as outfile:
            outfile.write(submission_xml.toprettyxml(indent='  ',
                                                     encoding='UTF-8'))

        self.submission_xml_fp = fp

    def write_all_xml_files(study_fp, sample_fp, experiment_fp, run_fp,
                            submission_fp, action):
        """Write all XML files needed for an EBI submission using current data

        Parameters
        ----------
        study_fp : str
            The filepath to which the study XML will be written
        sample_fp : str
            The filepath to which the sample XML will be written
        experiment_fp : str
            The filepath to which the experiment XML will be written
        run_fp : str
            The filepath to which the run XML will be written
        submission_fp : str
            The filepath to which the submission XML will be written
        action : {'ADD', 'VALIDATE', 'UPDATE'}
            What action to take when communicating with EBI

        Notes
        -----
        If any of the filepaths point to an existing file, it will be
        overwritten
        """
        self.write_study_xml(study_fp)
        self.write_sample_xml(sample_fp)
        self.write_experiment_xml(experiment_fp)
        self.write_run_xml(run_fp)
        self.write_submission_xml(submission_fp, action)

    def add_samples_from_templates(self, sample_template, prep_templates,
                                   per_sample_fastq_dir):
        """
        Parameters
        ----------
        sample_template : file
        prep_templates : list of file
        per_sample_fastq_dir : str
            Path to the direcotry containing per-sample FASTQ files containing
            The sequence labels should be:
            ``SampleID_SequenceNumber And Additional Notes if Applicable``
        """
        for sample in iter_file_via_list_of_dicts(sample_template):
            sample_name = sample.pop('sample_name')
            taxon_id = sample.pop('taxon_id', None)
            description = sample.pop('description', None)

            self.add_sample(sample_name, taxon_id=taxon_id,
                            description=description,
                            empty_value='unknown', **sample)

        for prep_template in prep_templates:
            for prep in iter_file_via_list_of_dicts(prep_template):
                sample_name = prep.pop('sample_name')
                platform = prep.pop('platform')
                experiment_design_description = prep.pop(
                    'experiment_design_description')
                library_construction_protocol = prep.pop(
                    'library_construction_protocol')

                file_path = join(per_sample_fastq_dir, sample_name+'.fastq')
                self.add_sample_prep(sample_name, platform, 'fastq',
                                     file_path, experiment_design_description,
                                     library_construction_protocol,
                                     empty_value='unknown', **prep)

    @classmethod
    def from_templates_and_demux_fastq(
            cls, study_id, study_title, study_abstract, investigation_type,
            sample_template, prep_templates, demux_seqs_fp, output_dir,
            **kwargs):
        """Generate an ``EBISubmission`` from templates and a sequence file

        Parameters
        ----------
        study_id : str
        study_title : str
        study_abstract : str
        investigation_type : str
        sample_template : file
        prep_templates : list of file
        demux_seqs_fp : str
            Path to FASTQ File containing the demultiplexed sequences in the
            format output by QIIME. Namely, the sequence labels should be:
            ``SampleID_SequenceNumber And Additional Notes if Applicable``
        output_dir : str
            The directory to output the per-sample FASTQ files. It will be
            created it if does not already exist. The files will be named
            <sample_name>.fastq
        """
        if not exists(output_dir):
            mkdir(output_dir)

        # generate the per-sample FASTQ files
        with open(demux_seqs_fp, 'U') as demux_seqs:
            split_helper(demux_seqs, output_dir, sequence_buffer_size=1000)

        # initialize the EBISubmission object
        submission = cls(study_id, study_title, study_abstract,
                         investigation_type, **kwargs)

        submission.add_samples_from_templates(
            sample_template, prep_templates, output_dir)

        return submission

    @classmethod
    def from_templates_and_per_sample_fastqs(cls, study_id, study_title,
                                             study_abstract,
                                             investigation_type,
                                             sample_template, prep_templates,
                                             per_sample_fastq_dir,
                                             **kwargs):
        """Generate an ``EBISubmission`` from templates and FASTQ files

        Parameters
        ----------
        study_id : str
        study_title : str
        study_abstract : str
        investigation_type : str
        sample_template : file
        prep_templates : list of file
        per_sample_fastq_dir : str
            Path to the direcotry containing per-sample FASTQ files containing
            The sequence labels should be:
            ``SampleID_SequenceNumber And Additional Notes if Applicable``

        Notes
        -----
        - kwargs will be passed directly to the ``EBISubmission`` constructor,
          which will add them as key-value pairs to the study attributes
          section of the submission
        """
        # initialize the EBISubmission object
        submission = cls(study_id, study_title, study_abstract,
                         investigation_type, **kwargs)

        submission.add_samples_from_templates(sample_template,
                                              prep_templates,
                                              per_sample_fastq_dir)

        return submission
