#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 16:57:41 2020

@author: bkim
"""

def S3act(ptn):
    
    orbit = []
    S3 = [[0,1,2], [0,2,1], [2,1,0], [1,0, 2], [1,2,0], [2,0,1]   ]
    for sigma in S3:
        temp = [ptn[sigma[0]], ptn[sigma[1]], ptn[sigma[2]]]
        if temp not in orbit:
            orbit.append(temp)
    return orbit

    

def tri_ptns(n):
    
    S = [2*i-1 for i in range(1, (n+1)//2 +1 ) ]
    PP = [ [ ] ] 
            
    while S:
        k = S.pop()
        PP2 = PP[:]
        for ptn in PP2:
            temp = ptn[:]+[k]
            temp.sort(reverse=True)
            if sum(temp) <=n:
                PP.append(temp)
    
    TP =[ [pi1, pi2, pi3] for pi1 in PP for pi2 in PP for pi3 in PP]
    TP =[ ptn for ptn in TP if sum(sum(pi) for pi in ptn) ==n ]
    
    TD = []
    while TP:
        ptn = TP.pop()
        orbit = S3act(ptn)
        if len(orbit)==3:
            ind =[i for i, pi in enumerate(orbit) if pi[1]==pi[2]]
            TD.append(orbit[ind[0]])
        else:
            TD.append(orbit[0])
        for ptn in orbit[1:]:
            if ptn in TP:
                TP.remove(ptn)
            
                
    return TD

class xi_map:
    
    def __init__(self, n):
        self.n = n
        self.domain = tri_ptns(n)
        self.codomain = tri_ptns(n+1)
        
    def order(self, pi, mu):
        larger = pi[:]
        smaller = mu[:]
        if pi==mu or mu==[]:
            pass
        else:
            if sum(larger)<sum(smaller):
                larger, smaller = mu[:], pi[:]
            elif sum(larger) == sum(smaller):
                for i in range(1, 1+min(len(larger), len(smaller))):
                    if sum(larger[:i]) < sum(smaller[:i]):
                        larger, smaller = mu[:], pi[:]
                        break
        
        return larger, smaller
    
    def triple_order(self, ptn):
        la0, la1, la2 = ptn[0][:], ptn[1][:], ptn[2][:]
        temp_lg, temp_sm = self.order(la0, la1)
        lg, _= self.order(temp_lg, la2)
        _, sm = self.order(temp_sm, la2)
        md = [la0, la1, la2]
        md.remove(lg)
        md.remove(sm)        
        return lg, md[0], sm
    
    def case_check(self, ptn):
        return len(S3act(ptn)), sum(sum(1 for pt in pi if pt==1) for pi in ptn)
    
    def map_from(self, ptn):
        j, k = self.case_check(ptn)
        pi, mu, la = ptn[0][:], ptn[1][:], ptn[2][:]
        pre = [pi, mu, la]
        
        
        if [j, k] == [6, 0]:
            pi, mu, la = self.triple_order([pi, mu, la])
            pi += [1]
            return [pi, mu, la]
        
        elif [j,k] ==[3,0] or [j,k] == [1,0]:
            return [pi+[1], mu, la]
        
        elif [j,k] == [6,1] or [j,k] == [3,1]:
            for a in pre:
                if 1 in a:
                    a.remove(1)
                else:
                    a += [1]
            return [pi, mu, la]
        
        elif [j, k] == [6,2]:
            
            pi, mu, la = self.triple_order([pi, mu, la])
            if 1 in pi:
                pi.remove(1)
                pi[0] +=2
            else:
                pi.append(1)
            return [pi, mu, la]
        
        elif [j, k] == [3,2]:
            return [pi+[1], mu, la]
        
        elif [j,k] == [6,3]:
            pi, mu, la = self.triple_order(pre)
            for a in [pi, mu, la]:
                a.remove(1)
            pi[0]+=4
            return [pi, mu, la]
        
        elif ([j,k] == [3,3] or [j, k]==[1,3]) and pi != [1]:
            for a in pre:
                a.remove(1)
            pi[0] +=4
            return [pi, mu, la]
        
        elif [j,k] == [3,3] and pi ==[1]:
            for a in pre:
                a.remove(1)
            mu[0]+=2
            la[0]+=2
            return [pi, mu, la]
        
    
    def detailed(self):
        filename = "map_for_n_" + str(self.n) + ".txt"
        with open(filename, "a") as f:
            print("The map xi on O_3/S_3(%d) to O_3/S_3(%d)" %(self.n, self.n+1), file=f)
            print("="*30, file = f)
            
            if self.well_defined_check():
                print("The images are in the codomain", file=f)
            
            if self.inj_check():
                print("The map xi is injective.", file =f)
                
            print("The detailed map is...", file = f)    
            for ptn in self.domain:
                j1, k1 = self.case_check(ptn)
                temp = self.map_from(ptn)
                j2, k2 = self.case_check(temp)
                
                print("""
    
    Partition %s 
    (in TD_{%d, %d})
    maps to 
    %s
    (in TD_{%d, %d})
                      """
                      %(ptn, j1, k1, temp, j2, k2)
                      , file =f)
                
            pass
        
    def in_check(self, ptn, ptn_set):
        temp = S3act(ptn)
        checker = False
        for a in temp:
            if a in ptn_set:
                checker = True
                re_ptn = a[:]
                break
            
        if checker:
            return re_ptn
        else:
            pass
        
        
    def well_defined_check(self):
        checker = True
        for ptn in self.domain:
            image = self.map_from(ptn)
            if self.in_check(image, self.codomain) ==None:
                print("something wrong")
                print(image)
                checker = False
        
        if checker:
            print("Images are in the codomain.")
        
        return checker
            
    
    def inj_check(self):
        
        range_set = []
        checker= True
        for ptn in self.domain:
            image = self.map_from(ptn)
            if self.in_check(image, range_set) == None:
                range_set.append(image)
            else:
                print("something wrong")
                print(ptn)
                print("maps to")
                print(image)
                checker = False
        if checker:
            print("It is injective.")
            
        return checker
                
                
            
                    
                


        
        
        
            
        
            
                



def main(n):
    xi = xi_map(n)
    xi.detailed()


if __name__=='__main__':
    for n in range(4, 31):
        main(n)
    