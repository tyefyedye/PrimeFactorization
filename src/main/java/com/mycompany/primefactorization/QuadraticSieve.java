package com.mycompany.primefactorization;

/**
 * Cesar Lobaton
 * CSCI 780 - Cryptography
 * RSA Project
 */

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;

public class QuadraticSieve extends PrimeFactoring {
    
    public QuadraticSieve(){ }
    
    public BigInteger getFactor(long b, BigInteger n){
        long bound; ArrayList<Long> primes;
        BigInteger p, pMinusOne, qr, a;
        ArrayList<BigInteger> factorBase = new ArrayList();
        ArrayList<Root> roots = new ArrayList();
        TreeMap<BigInteger, int[]> smoothNumbers = new TreeMap();
        
        System.out.println("Calculating bound...");
        if (b == -1) bound = getBound(n) + 1; else bound = b;
        System.out.println("Bound B = " + bound);
        System.out.println("Generating primes <= B...");
        primes = sieveOfEratosthenes(bound + 1);
        
        System.out.println("Reducing primes to quadratic residues over N only...");
        for (long prime : primes){
            p = BigInteger.valueOf(prime);
            pMinusOne = p.subtract(BI_ONE);
            qr = n.modPow(pMinusOne.divide(BI_TWO), p);
            if (qr.equals(BI_ONE)) factorBase.add(p);
        }
        int numFactors = factorBase.size();
        System.out.println(factorBase.size() + " quadratic residues found.");
        
        //a = n.sqrt().add(BI_ONE);
        a = n.sqrt();
        
        System.out.println("Getting roots over factor base...");
        BigInteger[] squareRoots = new BigInteger[2];
        BigInteger r1, r2;
        for (BigInteger fb : factorBase){
            if (fb.equals(BI_TWO)) roots.add(new Root(1, 2));
            else {
                squareRoots = tonelliShanks(n, fb);
                r1 = squareRoots[0].subtract(a).mod(fb);
                r2 = squareRoots[1].subtract(a).mod(fb);
                roots.add(new Root(r1.longValue(), fb.longValue()));
                roots.add(new Root(r2.longValue(), fb.longValue()));
            }
        }
        System.out.println("Roots found");
        
        System.out.println("Finding smooth numbers...");
        int numSmoothNumbers = 0, i = 0, j = 0;
        BigInteger root, f, t, r;
        int[] pExp;
        int[][] matrix = new int[numFactors][numFactors];
        while (numSmoothNumbers < numFactors){
            root = BigInteger.valueOf(roots.get(i % numFactors).getRoot());
            t = polynomial(a, root, n);
            r = t;
            if(!smoothNumbers.containsKey(t)){
                j = 0;
                pExp = new int[numFactors];
                Arrays.fill(pExp, 0);
                
                while (j < numFactors && !r.equals(BI_ONE)){
                    f = factorBase.get(j);
                    if (!r.mod(f).equals(BI_ZERO)) j++;
                    else { r = r.divide(f); pExp[j] = pExp[j] + 1; }
                }
                if (j < numFactors){
                    smoothNumbers.put(t, pExp);
                    matrix[numSmoothNumbers] = pExp;
                    numSmoothNumbers++;
                    if (numSmoothNumbers % 100 == 0)
                        System.out.println(numSmoothNumbers + " smooth numbers found");
                } else smoothNumbers.put(t, null);
            }
            roots.set(i % numFactors, roots.get(i % numFactors).incrementRoot());
            i++;
        }
        
        System.out.println(numSmoothNumbers + " total smooth numbers found");
        
        int[] arr;
        for(BigInteger bi : smoothNumbers.keySet()){
            arr = smoothNumbers.get(bi);
            if (arr != null) System.out.println(bi + ", " + Arrays.toString(arr));
        }
        
        return BI_ONE;
    }
    
    private BigInteger polynomial(BigInteger a, BigInteger i, BigInteger n) {
        BigInteger ai = a.add(i);
        return ai.pow(2).subtract(n);
    }
    
    
    public record Root (long root, long increment) {
        public Root incrementRoot() { return new Root(root + increment, increment); }
        
        public long getRoot() { return root; }
        
        public long getIncrement() { return increment; }
        
        @Override
        public String toString() { return "[" + root + ", " + increment + "]"; }
    }
}
