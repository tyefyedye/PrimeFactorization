package com.mycompany.primefactorization;

/**
 * Cesar Lobaton
 * CSCI 780 - Cryptography
 * RSA Project
 */

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;

public class QuadraticSieve extends PrimeFactoring {
    
    public QuadraticSieve(){ }
    
    public BigInteger getFactor(long b, BigInteger n){
        long bound; ArrayList<Long> primes;
        BigInteger p, pMinusOne, qr, a, t, r;
        BigInteger[] squareRoots = new BigInteger[2], pos = new BigInteger[2];
        ArrayList<BigInteger> factorBase = new ArrayList();
        TreeMap<BigInteger, ArrayList<Integer>> smoothNumbers = new TreeMap();
        
        if(b == -1) bound = getBound(n); else bound = b;
        primes = sieveOfEratosthenes(bound);
        
        for (long prime : primes){;
            p = BigInteger.valueOf(prime);
            pMinusOne = p.subtract(BI_ONE);
            qr = n.modPow(pMinusOne.divide(BI_TWO), p);
            if (qr.equals(BI_ONE)) factorBase.add(p);
        }
        int numFactors = factorBase.size();
        a = n.sqrt().add(BI_ONE);
        
        for (BigInteger fb : factorBase){
            System.out.println("Checking " + fb);
            if (fb.equals(BI_TWO)){
                pos[0] = BI_ZERO; pos[1] = BI_ONE;
            } else {
                squareRoots = tonelliShanks(n, fb);
                pos[0] = squareRoots[0].subtract(a).mod(fb);
                pos[1] = squareRoots[1].subtract(a).mod(fb);
            }

            int i = 0, count = 0;
            BigInteger j, pf, limit = new BigInteger("1000000");
            ArrayList<Integer> pExp;
            for(BigInteger ps : pos){
                count = 0; j = BI_ONE;
                while (count < 5 && !j.equals(limit)){
                    i = 0;
                    pExp = new ArrayList(Collections.nCopies(numFactors, 0));
                    t = polynomial(a.add(ps.multiply(j)), n);
                    r = t;
                    while(i < numFactors && !r.equals(BI_ONE)){
                        pf = factorBase.get(i);
                        if (!r.mod(pf).equals(BI_ZERO)) i++;
                        else {
                            r = r.divide(pf);
                            pExp.set(i, pExp.get(i) + 1);
                        }
                    }
                    if (i < numFactors && !smoothNumbers.containsKey(t)) {
                        smoothNumbers.put(t, pExp);
                    }
                    j = j.add(BI_ONE);
                }
            }
        }
        
        ArrayList<Integer> arr =  new ArrayList();
        for(BigInteger bi : smoothNumbers.keySet()){
            arr = smoothNumbers.get(bi);
            System.out.println(bi + ", " + arr);
        }
        
        return BI_ONE;
    }
    
    private BigInteger polynomial(BigInteger a, BigInteger n) {
        return a.pow(2).subtract(n);
    }
        
}
