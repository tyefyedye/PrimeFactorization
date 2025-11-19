package com.mycompany.primefactorization;

/**
 * Cesar Lobaton
 * CSCI 780 - Cryptography
 * RSA Project
 */

import org.jscience.mathematics.number.ModuloInteger;
import org.jscience.mathematics.number.LargeInteger;
import org.jscience.mathematics.vector.DenseVector;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;

public class QuadraticSieve extends PrimeFactoring {
    private ArrayList<Integer> dependents = new ArrayList();
    
    public QuadraticSieve(){ }
    
    public BigInteger getFactor(long b, BigInteger n){
        long bound; ArrayList<Long> primes;
        BigInteger p, pMinusOne, qr, a;
        ArrayList<BigInteger> factorBase = new ArrayList();
        ArrayList<Root> roots = new ArrayList();
        TreeMap<BigInteger, BigInteger> smoothNumbers = new TreeMap();
        TreeMap<BigInteger, double[]> smoothNumbersPrimeFactors = new TreeMap();
        
        System.out.println("Calculating bound...");
        if (b == -1) bound = getBound(n); else bound = b;
        System.out.println("Bound B = " + bound);
        System.out.println("Generating primes <= B...");
        primes = sieveOfEratosthenes(bound + 1);
        System.out.println(primes.size() + " primes found.");
        
        System.out.println("Reducing primes to quadratic residues over N only...");
        for (long prime : primes){
            p = BigInteger.valueOf(prime);
            pMinusOne = p.subtract(BI_ONE);
            qr = n.modPow(pMinusOne.divide(BI_TWO), p);
            if (qr.equals(BI_ONE)) factorBase.add(p);
        }
        int numFactors = factorBase.size();
        System.out.println(factorBase.size() + " quadratic residues found.");
        
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
        double[] pExp; 
        while (numSmoothNumbers < numFactors){
            root = BigInteger.valueOf(roots.get(i % numFactors).getRoot());
            t = polynomial(a, root, n);
            r = t;
            if(!smoothNumbers.containsKey(t)){
                j = 0;
                pExp = new double[numFactors];
                Arrays.fill(pExp, 0);
                
                while (j < numFactors && !r.equals(BI_ONE)){
                    f = factorBase.get(j);
                    if (!r.mod(f).equals(BI_ZERO)) j++;
                    else { r = r.divide(f); pExp[j] = (pExp[j] + 1); }
                }
                if (j < numFactors){
                    smoothNumbers.put(t, a.add(root));
                    smoothNumbersPrimeFactors.put(t, pExp);
                    numSmoothNumbers++;
                    if (numSmoothNumbers % 100 == 0)
                        System.out.println(numSmoothNumbers + " smooth numbers found");
                }
            }
            roots.set(i % numFactors, roots.get(i % numFactors).incrementRoot());
            i++;
        }
        System.out.println(numSmoothNumbers + " total smooth numbers found");
        
        System.out.println("Building matrix...");
        int k = 0; double[] arr;
        double[][] matrix = new double[numFactors][numFactors];
        for(BigInteger sn : smoothNumbers.keySet()){
            arr = smoothNumbersPrimeFactors.get(sn);
            if (arr != null) { matrix[k] = toModuloArray(arr, 2); k++; } 
        };
        
        System.out.println("Performing Gaussian Elimination on matrix...");
        matrix = gaussianElimination(matrix, numSmoothNumbers);
        System.out.println("Gaussian Elimination complete");
        
        return BI_ONE;
    }
    
    private BigInteger polynomial(BigInteger a, BigInteger i, BigInteger n) {
        BigInteger ai = a.add(i);
        return ai.pow(2).subtract(n);
    }
    
    private double[][] gaussianElimination(double[][] matrix, int length){
        ModuloInteger.setModulus(LargeInteger.valueOf(2));
        int ix; DenseVector<ModuloInteger> c, r, a;
        ArrayList<DenseVector<ModuloInteger>> columns = new ArrayList<DenseVector<ModuloInteger>>();
        boolean markedRows[] = new boolean[length]; Arrays.fill(markedRows, false);

        for (int i = 0; i < length; i++) {
            ArrayList<ModuloInteger> ml = new ArrayList<ModuloInteger>();
            for (int j = 0; j < length; j++) {
                ml.add(matrix[j][i] == 1 ? MI_1 : MI_0);
            }
            columns.add(DenseVector.valueOf(ml));
        }
        
//        for(DenseVector dv : columns) System.out.println(dv);
//        System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        
        for (int i = 0; i < length; i++){
            ix = 0; c = columns.get(i);
            while (ix < length && !c.get(ix).equals(MI_1)) ix++;
            if (ix < length) {
                markedRows[ix] = true;
//                System.out.println("Row " + ix + " is marked");
                for (int j = 0; j < columns.size(); j++){
                    r = columns.get(j);
                    if (!r.equals(c) && r.get(ix).equals(MI_1)){
                        a = r.plus(c);
                        columns.set(j, a);
                    }
                        
                }
            }
//            for(DenseVector dv : columns) System.out.println(dv);
//            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        }
        
        double[][] newMatrix = new double[length][length];
        for (int i = 0; i < length; i++) {
            DenseVector<ModuloInteger> dv = columns.get(i);
            for (int j = 0; j < length; j++) {
                newMatrix[j][i] = dv.get(j).equals(MI_1) ? 1 : 0;
            }
            if (!markedRows[i]) dependents.add(i);
        }
        
        return newMatrix;
    }
    
    
    public record Root (long root, long increment) {
        public Root incrementRoot() { return new Root(root + increment, increment); }
        
        public long getRoot() { return root; }
        
        public long getIncrement() { return increment; }
        
        @Override
        public String toString() { return "[" + root + ", " + increment + "]"; }
    }
}
