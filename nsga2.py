import numpy as np
from individual import Individual

class NSGA_2:
    def __init__(self,pop_size, sampling, crossover,
        mutation, eliminate_duplicates):
        self.pop_size = pop_size
        self.pm = 0.8
        return
    


    def set_up(self,problem,benchmark):
        self.lower = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0])
        self.upper = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
       3, 3, 3, 3])

        self.problem = problem
        self.benchmark = benchmark
        self.len_dna = problem.n_var
        self.n_objs = problem.n_obj

    def init_pop(self,pop_size):
        self.population = np.empty(pop_size,dtype=object)
        for idx in range(len(self.population)):
            individual = Individual(True,self.lower,self.upper,self.n_objs,None)
            self.population[idx] = individual
        return self.population
    def pop_cross(self,P):
        new_Q = []
        i = 0
        while i < len(P):
            q1,q2 = P[i].cross_over(P[i],P[i+1],False,self.lower,self.upper,None,self.n_objs)
            new_Q.append(q1)
            new_Q.append(q2)
            i += 2
        return new_Q
    
    def pop_mutation(self,Q):
        for q in Q:
            
            if np.random.rand() < self.pm:
                q.mutation(self.lower,self.upper)

        return Q
    
    def select(self,R):
        R = sorted(R)
        # print("sorted -------------------------------------------")
        # for item in R:

            # item.print_info()
        return R[:self.pop_size]
    
    def calculate_f(self,population):
        # print(len(population),"------------")
        self.F = np.zeros((len(population),self.n_objs))
        for idx in range(len(population)):

            dna = np.reshape(population[idx].dna,(1,-1))
            # print(dna)
            f = self.benchmark.evaluate(dna)
            # print(idx,f)
            population[idx].f = f[0]
    def fast_non_dominated_sort(self,population):
        n = np.zeros(len(population))  
        S = []  
        F = [[]]  
        for i in range(len(population)):
            Sp = []
            n[i] = 0
            for j in range(len(population)):
                if i == j:
                    continue
                if population[i].dominate(population[j]):
                    Sp.append(j)
                elif population[j].dominate(population[i]):
                    n[i] += 1
            S.append(Sp)
            if  n[i] == 0:
                population[i].rank = 0
                F[0].append(i)

        i = 0
        while len(F[i]) != 0:
            Q = []  
            for p in F[i]:
                for q in S[p]:
                    n[q] -= 1
                    if n[q] == 0:
                        population[q].rank = i + 1
                        Q.append(q)
            i += 1
            F.append(Q)
        F.pop()


        F_obj = []
        for i in range(len(F)):
            l = []

            for j in range(len(F[i])):
                l.append(population[F[i][j]])
            F_obj.append(l)

        return F_obj  
    
    def crowding_distance_pop(self,F):
        for rank in F:
            self.crowding_distance_assignment(rank)

        return
    
    def crowding_distance_assignment(self,rank):
        
        pop_fitness = np.zeros((len(rank),self.n_objs))

        for idx in range(len(rank)):
            rank[idx].crowding_dist = 0
            # print(rank[idx].f)
            pop_fitness[idx] = rank[idx].f

        sorted_idx = np.argsort(pop_fitness,axis=0)

        for i in range(self.n_objs):
            fi_sort_idx = sorted_idx[:,i]
            fi_sort_p = [rank[i] for i in fi_sort_idx]
            inf = (1e15)
            fi_sort_p[0].crowding_dist = inf
            fi_sort_p[-1].crowding_dist = inf
            for j in range(1,len(fi_sort_p)-1):
                fi_sort_p[j].crowding_dist += (fi_sort_p[j + 1].f[i] - fi_sort_p[j - 1].f[i])/(fi_sort_p[-1].f[i] - fi_sort_p[0].f[i])
            
        return rank


    def run(self,max_iteration):
        res = np.zeros((100,(self.len_dna)))
        P = self.init_pop(self.pop_size)
        self.calculate_f(P)
        F = self.fast_non_dominated_sort(P)
        self.crowding_distance_pop(F)
        R = P
        
        for idx in range(len(res)):
            res[idx] = R[idx].dna
        res = res.astype('int')
        hv = self.benchmark.calc_perf_indicator(res, 'hv')
        
        print(hv)
        print("------------------------------------------------------")
        l = []
        maxnum = 0
        for iter in range(100):
            Q = self.pop_cross(P)
            Q = self.pop_mutation(Q)
            R = np.concatenate((P,Q))
            
            # print(res)
            for individual in R:
                individual.rank = -1
                individual.crowding_dist = -1
                individual.f = np.zeros((self.n_objs))
            self.calculate_f(R)

            F = self.fast_non_dominated_sort(R)
            self.crowding_distance_pop(F)
            # for item in R:
            #     item.print_info()
            P = self.select(R)

            res = np.zeros((len(P),len(P[0].dna)))
            for idx in range(len(P)):
                res[idx] = P[idx].dna

            res = res.astype('int')
            # print(res)
            
            rank0 = []
            for i in range(len(P)):
                if P[i].rank == 0:
                    rank0.append(P[i].dna)
            print(len(rank0),"length of rank0")
            hv = self.benchmark.calc_perf_indicator(res, 'hv')
            rank0 = np.array(rank0)
            rank0 = rank0.astype('int')
            hv2 = self.benchmark.calc_perf_indicator(rank0,'hv')
            if hv > maxnum:
                maxnum = hv
            # print(hv,"hv of all")
            l.append(maxnum)

        import matplotlib.pyplot as plt
        x = np.arange(100)
        y = np.array(l)
        plt.plot(x,y)
        plt.xlabel('iteration') 
        plt.ylabel('hv')   

        plt.show()
        print(maxnum)
        return res
    
        # res = res.astype('int')
        
        # for i in res:
        #     print(i)
        # hv = self.banchmark.calc_perf_indicator(res, 'hv')
        # print(hv)
        # for i in range(max_iteration):
        #     print(i)
        #     P = self.select(R)

        #     res = np.zeros((len(P),len(P[0].dna)))
            # for idx in range(len(P)):
            #     res[idx] = P[idx].dna
            # res = res.astype('int')
            # hv = self.banchmark.calc_perf_indicator(res, 'hv')
            # print(hv)

        #     Q = self.pop_cross(P)
        #     Q = self.pop_mutation(Q)
        #     R = np.concatenate((P,Q))
        #     print('len R',len(R))
        #     self.calculate_f(R)
        #     F = self.fast_non_dominated_sort(R)
        #     self.crowding_distance_pop(F)
            
        
        # res = res.astype('int')
        # P = self.select(R)
        # res = np.zeros((len(P),len(P[0].dna)))
        # for idx in range(len(P)):
        #     res[idx] = P[idx].dna
        
        # res = res.astype('int')
        # return res




