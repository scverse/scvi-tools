class classproperty:
    """Read-only class property decorator.

    Source: https://stackoverflow.com/questions/5189699/how-to-make-a-class-property
    """

    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)
