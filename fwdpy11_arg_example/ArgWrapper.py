class ArgWrapper(object):
    __gc_interval = None

    def __init__(self, gc_interval):
        self.gc_interval = gc_interval

    def __call__(self, generation, nodes, edges):
        return False

    @property
    def gc_interval(self):
        return self.__gc_interval

    @gc_interval.setter
    def gc_interval(self, value):
        try:
            int(value)
        except:
            raise ValueError("GC interval must be an integer")
        if value <= 0:
            raise ValueError("GC interval must be and integer > 0")
        self.__gc_interval = int(value)
