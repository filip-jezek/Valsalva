from collections import defaultdict

class ModelicaClass:


    
    def __init__(self, name = 'Empty', start_val = None, path = None, isRoot = False, properties = None):
        self.name = name
        self.children = dict()
        self.start_val = start_val
        self.properties = properties
        if path is None:
            self.full_path = name 
        else: 
            self.full_path = path + '.'+ name
        if len(self.children) > 0 and self.isValueType:
            raise Exception("Value type " + name + " cannot have children!")
    
    @property
    def isValueType(self):
        return self.start_val is not None

    
    def findNode(self, full_path):

        if '.' in full_path:
            node, rest = full_path.split('.', 1)
        else:
            node = full_path

        mc = self.children.get(node)
        if mc is None:
            return None
        elif mc.isValueType:
            return mc
        else: 
            return mc.findNode(rest)

    def buildChildTree(self, inps):
        # ill = inp.split('\n')
        for line in inps:
            # line = line.split(',', 1)
            if len(line) > 1:
                l = line[0]
                p = line[1]
            else:
                l = line
                p = None

            self.__attach(l, self, properties = p)

    def printObjectTree(self, indent_level = 0) -> str:
        if self.isValueType:
            fixed = 'true'
            val = str(self.start_val)
            
            if self.properties is not None:
                if 'guess' in self.properties: 
                    fixed = 'false'
                elif 'bool' in self.properties:
                    val = 'true' if self.start_val == 1 else 'false'
                elif '__V_PV_init' in self.properties:
                    val = "%e + settings.V_PV_init" % self.start_val
                elif 'param' in self.properties:
                    return '%s = %s' % (self.name, val)
                elif 'discreteInit' in self.properties:
                    # same as stateInit
                    pass

            return '%s(start = %s, fixed = %s)' % (self.name, val, fixed)
        else :
            return '\n' + (' ' * indent_level) + self.name + '(' + self.printChildren(indent_level = indent_level+2) + ')'

    def printChildren(self, indent_level = 0) -> str:
        arr = (', ').join([self.children[c].printObjectTree(indent_level = indent_level) for c in self.children])
        return arr


    @staticmethod
    def __attach(branch, trunk, properties = None):
        '''
        Insert a branch of objects on its trunk.
        Thanks to https://stackoverflow.com/questions/8484943/construct-a-tree-from-list-os-file-paths-python-performance-dependent
        '''
        parts = branch.split('.', 1)
        if len(parts) == 1:  # branch is a file
            # add value type
            mc = ModelicaClass(parts[0], start_val=0, path=trunk.full_path, properties = properties)
            trunk.children[parts[0]] = mc
        else:
            node, others = parts
            if node not in trunk.children:
                # insert a non-value node
                mc = ModelicaClass(node, path=trunk.full_path)
                trunk.children[node] = mc
            ModelicaClass.__attach(others, trunk.children[node], properties=properties)

    @staticmethod
    def BuildObjectTree(lines, root = 'Root'):
        root_mc = ModelicaClass(root)
        root_mc.buildChildTree(lines)
        return root_mc


    @staticmethod
    def BuildObjectTreeFromFile(filename, root = 'Root', ignoreLines = 0):
        lines = open(filename, 'r').read().splitlines()
        return ModelicaClass.BuildObjectTree(lines[ignoreLines:], root)


