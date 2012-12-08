

class Peekorator(object):
    
    def __init__(self, iterator):
        self.empty = False
        self.peek = None
        self._iterator = iterator
        try:
            self.peek = iterator.next()
        except StopIteration:
            self.empty = True

    def __iter__(self):
        return self

    def next(self):        
        """
        Return the self.peek element, or raise StopIteration 
        if empty
        """
        if self.empty:
            raise StopIteration()
        to_return = self.peek
        try:
            self.peek = self._iterator.next()
        except StopIteration:
            self.peek = None
            self.empty = True
        return to_return


