"""
parse the conf and rule files
"""
from ConfigParser import SafeConfigParser

class LowerSectionConfigParser(SafeConfigParser):
    def _read(self, fp, fpname):
        SafeConfigParser._read(self, fp, fpname)
        new_sections = self._dict()
        for a_section in self._sections:
            tmp = self._sections[a_section]
            new_sections[a_section.lower()] = tmp
        self._sections = new_sections

class Section():
    def __init__(self, item_list):
        for i in item_list:
            self.__dict__[i] = None
    def __str__(self):
        return repr(self.__dict__)

class ConfBase():
    def __init__(self, conf):
        self.__conf__ = LowerSectionConfigParser()
        self.__conf__.read(conf)
        for section_name in dir(self):
            section = getattr(self, section_name)
            if isinstance(section, Section):
                for option_name in section.__dict__:
                    section.__dict__[option_name] = self.__conf__.get(section_name, option_name)

class RuleBase():
    def __init__(self, conf, dataid, controlnumber, treatnumber):
        self.__conf__ = LowerSectionConfigParser()
        self.__conf__.read(conf)
        # use default node to set datasetid
        self.__conf__.set('DEFAULT',
                          'DatasetID', dataid)
        # Convert treat, control in to list
        # get replicates number
        for section_name in dir(self):
            get_raw = lambda opt: self.__conf__.get(section_name, opt, 1)
            get_expand = lambda opt: self.__conf__.get(section_name, opt, 0)
            def get_fmt(opt):
                get_rep_expand = lambda raw_str, rep_cnt: map(lambda x:self.__conf__.get(section_name, opt, 0, {raw_str: str(x+1)}),
                                                         range(rep_cnt))
                if '(treat_rep)' in get_raw(opt):
                    return get_rep_expand("treat_rep", treatnumber)
                elif '(control_rep)' in get_raw(opt):
                    return get_rep_expand("control_rep", controlnumber)
                else:
                    return get_expand(opt)
            section = getattr(self, section_name)
            if isinstance(section, Section):
                for option_name in section.__dict__:
                    section.__dict__[option_name] = get_fmt(option_name)
