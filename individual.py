import numpy as np
class Individual:
    def __init__(self, random_init, lower, upper, dna,n_objs):
        if random_init:
            self.dna = np.random.randint(lower, upper)
        else:
            self.dna = dna

        self.rank = -1
        self.crowding_dist = -1
        self.f = np.zeros((n_objs))

    def print_info(self):
        print("dna",self.dna," rank",self.rank,"F",self.f," crowding distance",self.crowding_dist)

    def __lt__(self,other):
        if self.rank < other.rank:
            return True
        if self.rank == other.rank and self.crowding_dist > other.crowding_dist:
            return True
        else:
            return False
    
    def dominate(self,other):
        for idx in range(len(self.f)):
            if self.f[idx] > other.f[idx]:
                return False
        return True

    def mutation(self,lower,upper):
        idx = np.random.randint(0,len(self.dna))
        l = lower[idx]
        r = upper[idx]
        self.dna[idx] = np.random.rand()*(r-l)+r
        return self.dna
    
    def cross_over(self,p1,p2,random_init,lower,upper,dna,n_objs):
        p1_dna = p1.dna.copy()
        p2_dna = p2.dna.copy()

        idx = np.random.randint(1,len(p1_dna)-1)

        temp = p1_dna[idx:].copy()
        p1_dna[idx:] = p2_dna[idx:]
        p2_dna[idx:] = temp

        p1_new = Individual(random_init, lower,upper,p1_dna,n_objs)
        p2_new = Individual(random_init, lower,upper,p2_dna,n_objs)
        return p1_new, p2_new