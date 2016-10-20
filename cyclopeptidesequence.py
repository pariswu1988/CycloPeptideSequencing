peptide_masses =    {"A":  71,
                    "C": 103,
                    "D": 115,
                    "E": 129,
                    "F": 147,
                    "G":  57,
                    "H": 137,
                    "I": 113,
                    "K": 128,
                    "L": 113,
                    "M": 131,
                    "N": 114,
                    "P":  97,
                    "Q": 128,
                    "R": 156,
                    "S":  87,
                    "T": 101,
                    "V":  99,
                    "W": 186,
                    "Y": 163}

default_es = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def CyclicSpectrum(peptide):
    prefix_mass = [0]
    for i in range(1, len(peptide) + 1):
        temp = prefix_mass[i - 1] + peptide_masses[peptide[i - 1]]
        prefix_mass.append(temp)
    peptide_mass = prefix_mass[len(peptide)]
    cyclo_spectrum = [0]
    for i in range(0, len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            cyclo_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                cur = peptide_mass - (prefix_mass[j] - prefix_mass[i])
                cyclo_spectrum.append(cur)
    return sorted(cyclo_spectrum)

def LinearSpectrum(peptide):
    prefix_mass = [0]
    for i in range(1, len(peptide) + 1):
        temp = prefix_mass[i - 1] + peptide_masses[peptide[i - 1]]
        prefix_mass.append(temp)
    linear_spectrum = [0]
    for i in range(0, len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)

def TestCyclicSpectrum(peptide):
    a = []
    for pep in peptide:
        a.append(peptide_masses[pep])
    return a

def Expand(peptides, expansion_set = default_es):
    expanded_peptides = []
    for peptide in peptides:
        for key in expansion_set:
            new_peptide = peptide + key
            expanded_peptides.append(new_peptide)
    return expanded_peptides

def Mass(peptide):
    m = 0
    for pep in peptide:
        m += peptide_masses[pep]
    return m

def ParentMass(spectrum):
    return spectrum[-1]

def Inconsistent(peptide, spectrum):
    a = LinearSpectrum(peptide)
    s = spectrum[:]
    for pep in a[:]:
        if pep not in s:
            return True
        s.remove(pep)
    return False

def CyclopeptidesSequencing(spectrum):
    peptides = [""]
    ans = []
    Testpeptides = Expand(peptides)
    es = []
    for p in Testpeptides[:]:
        if not Inconsistent(p, spectrum):
            es.append(p)
    while peptides:
        peptides = Expand(peptides, es)
        for p in peptides[:]:
            if Mass(p) == ParentMass(spectrum):
                if CyclicSpectrum(p) == spectrum:
                    temp = Spectrify(p)
                    ans.append(temp)
                peptides.remove(p)
            elif Inconsistent(p, spectrum):
                peptides.remove(p)
    ans = RemoveDupes(ans)
    return ans

def Spectrify(p):
    a = []
    for c in p:
        t = peptide_masses[c]
        a.append(t)
    return a

def RemoveDupes(l):
    s = []
    for i in l:
        if i not in s:
            s.append(i)
    return s

def PrettyPrint(specs):
    for s in specs:
        a = ""
        for n in s:
            a += "-"
            a += str(n)
        print(a[1:])


with open('cyclo_example.txt') as f:
    t = f.read().strip()
    a = []
    for n in t.split(" "):
        a.append(int(n))

PrettyPrint(CyclopeptidesSequencing(a))

