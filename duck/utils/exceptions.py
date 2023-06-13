class EquilibrationError (Exception):
        def __init__(self, msg="System is not equilibrated.", *args, **kwargs):
            super().__init__(msg, *args, **kwargs)