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
    
def logicOutput(P,Q,R):
    return (P>=(~Q|R))^(R|Q)

num = 3 #number of variables
varspace = 3
outspace = 8
print("P".rjust(varspace," "), "Q".rjust(varspace," "), "R".rjust(varspace," "), "output".rjust(outspace," "),end ="\n\n")
for i in range(1<<num):#2^number of variables
    varVal = []
    for x in range(num):
        result = logic()
        result.bool = bool(i&(1<<x))
        varVal.append(result)
    
    [print("T".rjust(varspace," "), end = " ") if x.bool \
     else print("F".rjust(varspace," "),end=" ") for x in varVal]
    if logicOutput(*varVal).bool:
        out = "T"
    else: out = "F"
    print(out.rjust(outspace," "))
    