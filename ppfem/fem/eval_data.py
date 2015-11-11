

class EvalData(object):
    def __init__(self, attr_dict):
        self.add_attributes(attr_dict)

    def add_attributes(self, attr_dict):
        for k in attr_dict.keys():
            self.add_attribute(k, attr_dict[k])

    def add_attribute(self, name, value):
        self.__setattr__(name, value)

    def remove_attribute(self, name):
        self.__delattr__(name)

    def remove_attributes(self, names):
        for name in names:
            self.remove_attribute(name)
