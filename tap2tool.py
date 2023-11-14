import binascii
import errno
import functools
import hashlib
import io
import signal
import sys
import os
import time

import chardet
import json
import numpy
import pyvo
from pyvo import registry as pyvo_registry
from pyvo import DALQueryError

import urllib
from urllib import request, parse

from astropy.io import fits

from astropy import units as u
from astropy.coordinates import SkyCoord


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('execution time: %r  %2.2f ms' % \
                  (method.__name__, (te - ts) * 1000))
        return result
    return timed


class NonIVOATapArchiveWrapper:

    def __init__(self, access_url, res_title="", short_name=""):
        self.access_url = access_url
        self.res_title = res_title
        self.short_name = short_name


class TapArchive:

    OBSCORE_TABLE = "ivoa.obscore"

    def __init__(self, archive_object):
        self.archive_object = archive_object
        self.initialized = False
        self.service = None
        self.tables = None

        self.initialized = self._initialize()

        if self.initialized:
            self._set_name()
            self._set_hash()

    def _initialize(self):
        try:
            self.access_url = self.archive_object.access_url
        except Exception as e:
            return False

        if self.access_url == "" or self.access_url is None:
            return False

        self._set_service()

        try:
            self._set_tables()
        except Exception as e:
            return False

        return True

    def _set_name(self):

        self.name = ''

        if self.archive_object.res_title != '' and self.archive_object.res_title is not None:
            self.name = self.archive_object.res_title
        elif self.archive_object.short_name != '' and self.archive_object.short_name is not None:
            self.name = self.archive_object.short_name
        else:
            self.name = self.access_url

    def _set_hash(self):
        if not self.name:
            self._set_name()

        self.hash = hashlib.md5(self.name.encode("utf-8")).hexdigest()

    def _set_service(self):
        self.service = pyvo.dal.TAPService(self.access_url)

    @timeout(5)
    def _set_tables(self):
        tables = []

        try:
            for table in self.service.tables:
                fields = []

                for table_field in table.columns:
                    field = TapArchive.Field(table_field.name,
                                             table_field.datatype.content,
                                             table_field.unit,
                                             table_field.description)

                    fields.append(field)

                archive_table = TapArchive.Table(table.name, table.type, fields)

                tables.append(archive_table)
        except Exception as e:
            print(e)

        self.tables = tables

    def get_name(self):
        if not self.name:
            self._set_name()

        return self.name

    def get_hash(self):
        if not self.hash:
            self._set_hash()

        return self.hash

    def _is_ivoa_compliant(self):
        is_compliant = False

        try:
            for table in self.service.tables:
                if str(table.name).lower() == TapArchive.OBSCORE_TABLE:
                    is_compliant = True
                    return is_compliant
        except Exception as e:
            print(e)

        return is_compliant

    def is_ivoa(self):
        return self._is_ivoa_compliant()

    class Table:

        def __init__(self, name, type, fields):
            self.name = name
            self.type = type
            self.fields = fields

    class Field:

        def __init__(self, name, type, unit, description):
            self.name = name
            self.type = type
            self.unit = unit
            self.description = description


class QueryBuildersBlock:

    BLOCK_NAME = "query_type"
    BLOCK_LABEL = "ADQL Query Selection"

    TABLE_SELECTION_BLOCK_NAME = "table_selection"
    TABLE_SELECTION_BLOCK_LABEL = "Table selection"

    GENERIC_QUERY_NAME = "none"
    GENERIC_QUERY_LABEL = "No specific query (first n files from archive)"

    IVOA_BUILDER_NAME = "obscore_query"
    IVOA_BUILDER_LABEL = "IVOA obscore table query builder"

    RAW_QUERY_NAME = "raw_query"
    RAW_QUERY_LABEL = "Raw ADQL query builder"

    ACCESS_URL_INPUT_NAME = "access_url"
    ACCESS_URL_INPUT_LABEL = "Url download field"
    ACCESS_URL_INPUT_HELP = "Field containing the download url of the file"

    def __init__(self, archive_list):
        self.archive_list = archive_list
        self.builder_selection_block = ""
        self.builders_block = ""

    def generate_block(self):
        self.builder_selection_block = self.generate_builder_selection_block()
        self.builders_block = self.generate_builders_block()

    def generate_builder_selection_block(self):
        builder_selection_block = ""
        options = []

        option = ToolOption(QueryBuildersBlock.GENERIC_QUERY_NAME, QueryBuildersBlock.GENERIC_QUERY_LABEL)
        options.append(option)

        option = ToolOption(QueryBuildersBlock.IVOA_BUILDER_NAME, QueryBuildersBlock.IVOA_BUILDER_LABEL)
        options.append(option)

        option = ToolOption(QueryBuildersBlock.RAW_QUERY_NAME, QueryBuildersBlock.RAW_QUERY_LABEL)
        options.append(option)

        for archive in self.archive_list:
            option = ToolOption(archive.get_hash(), archive.get_name())
            options.append(option)

        builder_selection_block = ToolSelect(QueryBuildersBlock.BLOCK_NAME, QueryBuildersBlock.BLOCK_LABEL, options)

        return builder_selection_block

    def generate_builders_block(self):
        builders_block = ""

        for archive in self.archive_list:

            xml_params = ""

            query_builder_table_selection_block = self.generate_builder_table_selection_block(archive.tables)
            query_builder_params = self.generate_builder_params_block(archive.tables)

            xml_params = query_builder_table_selection_block.get_xml()

            for qbp in query_builder_params:
                xml_params += qbp.get_xml()

            query_builder_block = ToolWhenNested(archive.get_hash(), xml_params)

            builders_block += query_builder_block.get_xml()

        return builders_block

    def generate_builder_params_block(self, archive_tables):
        query_builders_params = []

        for table in archive_tables:
            fields_params = []

            download_url_param_options = []

            download_url_param = ToolParam(QueryBuildersBlock.ACCESS_URL_INPUT_NAME,
                                           QueryBuildersBlock.ACCESS_URL_INPUT_LABEL,
                                           "text",
                                           QueryBuildersBlock.ACCESS_URL_INPUT_HELP)

            fields_params.append(download_url_param)

            for field in table.fields:
                download_url_param_option = ToolOption(field.name, field.name)
                download_url_param_options.append(download_url_param_option)

                field_param = ToolParam(field.name, field.name, data_type=field.type, help=field.description)
                fields_params.append(field_param)

            tool_select_name = "download_url"
            download_url_select = ToolSelect(tool_select_name,
                                             "Download field",
                                             download_url_param_options)

            fields_params.insert(0, download_url_select)

            table_params = ToolWhen(table.name, fields_params)

            query_builders_params.append(table_params)

        return query_builders_params

    def generate_builder_table_selection_block(self, archive_tables):
        builder_table_selection_block = ""
        options = []

        for table in archive_tables:
            option = ToolOption(table.name, table.name)
            options.append(option)

        builder_table_selection_block = ToolSelect(
            QueryBuildersBlock.TABLE_SELECTION_BLOCK_NAME,
            QueryBuildersBlock.TABLE_SELECTION_BLOCK_LABEL,
            options)

        return builder_table_selection_block

    def get_xml(self):
        return self.builder_selection_block.get_xml() + self.builders_block


class ToolParam:

    DEFAULT_TYPE = "text"

    allowed_type_list = [
        "int",
        "text",
        "float",
        "select"
    ]

    def __init__(self, name, label, data_type=None, help=None):
        self.name = name
        self.label = label

        if data_type:
            self.data_type = data_type
            if not ToolParam._is_type_valid(data_type):
                self.data_type = ToolParam.DEFAULT_TYPE
        else:
            self.data_type = ToolParam.DEFAULT_TYPE

        self.help = help

    @staticmethod
    def _is_type_valid(data_type):
        if data_type in ToolParam.allowed_type_list:
            return True
        else:
            return False

    def get_xml(self):
        if self.data_type and self.help:
            return f'<param name="{self.name}" type="{self.data_type}" label="{self.label}" help="{self.help}" />'
        elif self.data_type:
            return f'<param name="{self.name}" type="{self.data_type}" label="{self.label}" />'
        elif self.help:
            return f'<param name="{self.name}" type="{ToolParam.DEFAULT_TYPE}" label="{self.label}" help="{self.help}" />'
        else:
            return f'<param name="{self.name}" type="{ToolParam.DEFAULT_TYPE}" label="{self.label}" />'


class ToolOption:

    def __init__(self, value, name):
        self.value = value
        self.name = name

    def get_xml(self):
        return f'<option value="{self.value}">{self.name}</option>'


class ToolSelect:

    def __init__(self, name, label, options):
        self.name = name
        self.label = label
        self.options = options

    def get_xml(self):
        select = ""

        select_top = f'<param name="{self.name}" type="select" label="{self.label}">'
        select_bot = '</param>'

        select += select_top + '\n'

        for option in self.options:
            select += '\t' + option.get_xml() + '\n'

        select += select_bot + '\n'

        return select


class ToolWhen:

    def __init__(self, value, params):
        self.value = value
        self.params = params

    def get_xml(self):
        when_statement = ""

        when_top = f'<when value="{self.value}">' + '\n'
        when_bot = '</when>'

        for param in self.params:
            when_top += '\t' + param.get_xml() + '\n'

        when_statement = when_top + when_bot + '\n'

        return when_statement


class ToolWhenNested:

    def __init__(self, value, xml_params):
        self.value = value
        self.xml_params = xml_params

    def get_xml(self):
        when_statement = ""

        when_top = f'<when value="{self.value}">' + '\n'
        when_top += ToolConditional.CONDITIONAL_TAG_OPEN

        when_bot = ToolConditional.CONDITIONAL_TAG_CLOSE
        when_bot += '</when>'

        when_top += self.xml_params + '\n'

        when_statement = when_top + when_bot + '\n'

        return when_statement


class ToolConditional:

    CONDITIONAL_TAG_OPEN = "<conditional> \n"
    CONDITIONAL_TAG_CLOSE = "</conditional> \n"

    def __init__(self, conditional_name, xml_block):
        self.conditional_name = conditional_name
        self.xml_block = xml_block

    def get_conditional_block(self):
        conditional_block = ToolConditional.CONDITIONAL_TAG_OPEN \
                            + self.xml_block \
                            + ToolConditional.CONDITIONAL_TAG_CLOSE

        return conditional_block


class FileHandler:

    def __init__(self):
        pass

    @staticmethod
    def write_file_to_output(file, output, write_type="w"):
        with open(output, write_type) as file_output:
            file_output.write(file)


@timeout(5)
def get_archive_tables(archive):
    return archive.search("SELECT TOP 100 * FROM ivoa.obscore")


@timeout(5)
def get_archive_fields(archive):
    return archive.search("SELECT TOP 100 * FROM ivoa.obscore")


@timeout(5)
def _set_archive_tables(archive_service):

    tables = []

    for table in archive_service.tables:
        archive_table = {
            'name': table.name,
            'type': table.type,
            'fields': None
        }

        fields = []

        for table_field in table.columns:
            field = {
                'name': table_field.name,
                'description': table_field.description,
                'unit': table_field.unit,
                'datatype': table_field.datatype.content
            }

            fields.append(field)

        archive_table['fields'] = fields

        tables.append(archive_table)

        print(tables)


@timeout(50)
def _set_archive_tables_summary(archive_service):

    tables = []

    print(archive_service.tables)

    try:
        for table in archive_service.tables:
            print(table.name)

            for field in table.columns:
                print(field.name)
    except Exception as e:
        print(e)

    #print(tables)


def _create_xml_conditional_archive(archive_name, archive_hash):
    option = f'<option value="{archive_hash}">{archive_name} query builder</option>'
    print(option)
    return option


def _create_xml_query_builders_list(conditional_archive_list):
    builders_list = ""

    builders_list_top = '<param name="query_type" type="select" label="ADQL Query Selection">'

    builders_list_none = '<option value="none">No specific query (first n files from archive)</option>'
    builders_list_ivoa = '<option value="obscore_query">IVOA obscore table query builder</option>'
    builders_list_raw = '<option value="raw_query">Raw ADQL query builder</option>'

    builders_list_bot = '</param>'

    builders_list += builders_list_top + '\n'
    builders_list += '\t' + builders_list_none + '\n'
    builders_list += '\t' + builders_list_ivoa + '\n'
    builders_list += '\t' + builders_list_raw + '\n'

    for conditional_archive in conditional_archive_list:
        builders_list += '\t' + conditional_archive + '\n'

    builders_list += builders_list_bot + '\n'


def _create_xml_archive_query_builder(archive_hash, archive_tables):
    query_builder = ""

    query_builder_top = '<when value = "coordinates">'

    for archive_table in archive_tables: 
    
        param = '<param name = "ra" type = "text" label = "Right ascension" />'

    query_builder_bot = '</when>'


def _get_archive_name(archive_detail):

    name = ''

    if archive_detail.res_title != '' and archive_detail.res_title is not None:
        name = archive_detail.res_title
    elif archive_detail.short_name != '' and archive_detail.short_name is not None:
        name = archive_detail.short_name
    else:
        name = archive_detail.access_url

    print(name)
    print(hashlib.md5(name.encode("utf-8")).hexdigest())

    return hashlib.md5(name.encode("utf-8")).hexdigest(), name


def get_non_ivoa_archives(archive_list):

    non_ivoa_archives = []

    #for i in range(archive_list.__len__()):
    for i in range(10):
        archive_record = archive_list.getrecord(i)
        archive = TapArchive(archive_record)
        if archive.initialized:
            if not archive.is_ivoa():
                non_ivoa_archives.append(archive)

    return non_ivoa_archives


def get_tap_archives():
    archive_list = pyvo.registry.search(servicetype="tap")

    return archive_list


def create_query_builders_xml_block():
    xml_query_builder_block = ""

    archive_list = get_tap_archives()
    curated_archive_list = get_non_ivoa_archives(archive_list)

    qbb = QueryBuildersBlock(curated_archive_list)

    return xml_query_builder_block


def create_xml_query_builders_block_file(file_path):
    file_content = ""
    file_content = create_query_builders_xml_block()

    FileHandler.write_file_to_output(file_content, file_path, "w")


def create_archive_object_non_ivoa(access_url):
    archive_object = NonIVOATapArchiveWrapper(access_url)

    return archive_object


def test_custom():
    archive_list = pyvo.registry.search(servicetype="tap")

    conditional_archive_list = []

    service = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")

    print(service)

    _set_archive_tables_summary(service)

    data = service.search("SELECT TOP 10 * FROM allwise.source")

    print(data.getrecord(1))

    # for i in range(archive_list.__len__()):
    for i in range(10):

        try:
            print(archive_list.getrecord(i).access_url)
            service = pyvo.dal.TAPService(archive_list.getrecord(i).access_url)
            _set_archive_tables_summary(service)
            archive_hash, archive_name = _get_archive_name(archive_list.getrecord(i))
            _create_xml_conditional_archive(archive_name, archive_hash)
            conditional_archive_list.append(_create_xml_conditional_archive(archive_name, archive_hash))
        except Exception as e:
            print(e)

    for conditional_archive in conditional_archive_list:
        print(conditional_archive)


def test_main():
    curated_tap_archives = []

    archives = get_tap_archives()

    curated_tap_archives = get_non_ivoa_archives(archives)

    #custom_archive = create_archive_object_non_ivoa("https://datalab.noirlab.edu/tap")
    #custom_archive = TapArchive(custom_archive)

    #curated_tap_archives.append(custom_archive)

    qbb = QueryBuildersBlock(curated_tap_archives)

    qbb.generate_block()

    builders_block_xml = qbb.get_xml()

    print(builders_block_xml)

    file_location = "/tool_build.xml"
  
    FileHandler.write_file_to_output(builders_block_xml, file_location, "w")


if __name__ == '__main__':

    test_main()
