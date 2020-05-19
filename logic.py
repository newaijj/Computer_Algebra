class logic:
    def __init__(self, b = None):
        self.bool = b
    def __xor__(self, other):
        result = logic()
        result.bool = self.bool and other.bool
        return result
    def __or__(self, other):
        result = logic()
        result.bool = self.bool or other.bool
        return result
    def __invert__(self):
        result = logic()
        result.bool = not self.bool
        return result
    def __ge__(self, other):
        result = logic()
        result.bool =(not self.bool) or other.bool
        return result
    def __str__(self):
        return str(self.bool)
    
P = logic(True)
Q = logic(True)
R = logic(False)
S = logic(False)

print("P^Q", P^Q)
print("P^R", P^R)
print("~P=>R", ~P>=R)
print("R=>P", R>=P)