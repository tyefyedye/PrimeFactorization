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
    private final BigInteger BI_THREE = BigInteger.valueOf(3);
    private final BigInteger BI_NINE = BigInteger.valueOf(9);
    private final BigInteger BI_TEN = BigInteger.valueOf(10);
    
    private TreeMap<BigInteger, BigInteger> smoothNumbers = new TreeMap();
    private TreeMap<BigInteger, double[]> smoothNumbersPrimeFactors = new TreeMap();
    private ArrayList<Integer> dependents = new ArrayList();
    private long lastBoundUsed;
    
    public QuadraticSieve(){ }
    
    public BigInteger getFactor(BigInteger n, long b, int range, int threshold, int minCandidates){
        if(n.equals(BI_NINE)) return BI_THREE;
        else if (n.compareTo(BI_TEN) == -1) return BI_TWO;
        
        long bound; ArrayList<Long> primes;
        BigInteger p, pMinusOne, qr, factor = BI_ONE;
        int totalSmoothNumbers;
        ArrayList<BigInteger> factorBase = new ArrayList();
        
        System.out.println("Calculating bound...");
        if (b == -1) bound = getBound(n); else bound = b;
        lastBoundUsed = bound;
        System.out.println("Bound B = " + bound);
        
        int numFactors = 0;
        while (numFactors < 7){
             System.out.println("Generating primes <= B...");
            primes = sieveOfEratosthenes(bound + 1);
            System.out.println(primes.size() + " primes found.");

            System.out.println("Reducing primes to quadratic residues over n...");
            for (long prime : primes){
                p = BigInteger.valueOf(prime);
                pMinusOne = p.subtract(BI_ONE);
                qr = n.modPow(pMinusOne.divide(BI_TWO), p);
                if (qr.equals(BI_ONE)) factorBase.add(p);
            }
            numFactors = factorBase.size();
            System.out.println(numFactors + " quadratic residues found.");
            if (numFactors < 7) {
                System.out.println("Not enough quadratic residues found.");
                factorBase.clear(); bound = bound + 10;
                System.out.println("Trying bound B = " + bound);
            }
        }
        
        System.out.println("Finding smooth numbers...");
        smoothNumbers.clear(); smoothNumbersPrimeFactors.clear();
        totalSmoothNumbers = getSmoothNumbers(n, factorBase, range, threshold, minCandidates);
        System.out.println(totalSmoothNumbers + " smooth numbers found");
       
        System.out.println("Building matrix...");
        int k = 0; double[] arr;
        double[][] matrix = new double[totalSmoothNumbers][numFactors];
        for(BigInteger sn : smoothNumbers.keySet()){
            arr = smoothNumbersPrimeFactors.get(sn);
            matrix[k] = arr; k++;
        }
        System.out.println("Matrix built");

        System.out.println("Performing Gaussian Elimination on matrix...");
        dependents.clear();
        matrix = gaussianElimination(matrix, totalSmoothNumbers, numFactors);
        System.out.println("Gaussian Elimination complete");

        System.out.println("Finding factor...");
        factor = findFactor(matrix, totalSmoothNumbers, numFactors, smoothNumbers, n);
        if(!factor.equals(BI_ONE) && !factor.equals(n)){
            System.out.println("Factor found!");
        }
        
        return factor;
    }
    
    public long incrementBound() { return lastBoundUsed + 10; }
    
    private int getSmoothNumbers(BigInteger n, ArrayList<BigInteger> factorBase,
            int range, int threshold, int minCandidates){
        int numFactors = factorBase.size(), usedRange = range, numCandidates = 999999;
        BigInteger a = n.sqrt(), ix, num, f;
        BigInteger[] polyNums = null; int[] sieve = null;
        ArrayList<BigInteger> candidates = new ArrayList();
        
        while (numCandidates > minCandidates){
            polyNums = new BigInteger[usedRange+1];
            sieve = new int[usedRange+1];
            for (int i = 0; i <= usedRange; i++){
                ix = BigInteger.valueOf(i);
                num = a.add(ix).modPow(BI_TWO, n);
                polyNums[i] = num;
                sieve[i] = num.bitLength();
            }

            int idxStart = factorBase.getFirst().equals(BI_TWO) ? 1 : 0;
            BigInteger[] squareRoots; int logF, idx;
            for (int i = idxStart; i < numFactors; i++){
                f = factorBase.get(i);
                squareRoots = tonelliShanks(n, f);
                logF = f.bitLength();
                for (BigInteger root : squareRoots){
                    idx = root.subtract(a).mod(f).intValue();
                    while (idx < usedRange) {
                        sieve[idx] = sieve[idx] - logF;
                        idx += f.intValue();
                    }
                }
            }

            for (int i = 0; i <= usedRange; i++)
                if(sieve[i] < threshold) candidates.add(polyNums[i]);
            
            numCandidates = candidates.size();
            if(numCandidates > 10000){
                System.out.println("Too many candidates! (" + numCandidates + ") Reducing range.");
                candidates = new ArrayList();
                usedRange = usedRange / 10;
            } else System.out.println(candidates.size() + " candidates found");
        }

        double[] pExp; int j, numPosSmoothNumbers = 0; BigInteger r, i;
        for(BigInteger c : candidates){
            if (!smoothNumbers.containsKey(c)){
                pExp = new double[numFactors];
                Arrays.fill(pExp, 0);
                j = 0;
                r = c;
                while (j < numFactors && !r.equals(BI_ZERO) && !r.equals(BI_ONE)){
                    f = factorBase.get(j);
                    if (!r.mod(f).equals(BI_ZERO)) j++;
                    else { r = r.divide(f); pExp[j] = (pExp[j] + 1) % 2; }
                }
                if (j < numFactors){
                    i = BigInteger.valueOf(getIndex(polyNums, c));
                    smoothNumbers.put(c, a.add(i));
                    smoothNumbersPrimeFactors.put(c, pExp);
                    numPosSmoothNumbers++;
                }
            }
        }
        return numPosSmoothNumbers;
    }
    
    private int getIndex(BigInteger[] nums, BigInteger numToSearch){
        for (int i = 0; i < nums.length; i++)
            if (nums[i].equals(numToSearch)) return i;
        return -1;
    }
    
    private double[][] gaussianElimination(double[][] matrix, int sm, int f){
        ModuloInteger.setModulus(LargeInteger.valueOf(2));
        int ix; DenseVector<ModuloInteger> c, r, a;
        ArrayList<DenseVector<ModuloInteger>> columns = new ArrayList();
        boolean markedRows[] = new boolean[sm]; Arrays.fill(markedRows, false);

        for (int i = 0; i < f; i++) {
            ArrayList<ModuloInteger> ml = new ArrayList();
            for (int j = 0; j < sm; j++) {
                ml.add(matrix[j][i] == 1 ? MI_1 : MI_0);
            }
            columns.add(DenseVector.valueOf(ml));
        }

        for (int i = 0; i < f; i++){
            ix = 0; c = columns.get(i);
            while (ix < sm && !c.get(ix).equals(MI_1)) ix++;
            if (ix < sm) {
                markedRows[ix] = true;
                for (int j = 0; j < f; j++){
                    r = columns.get(j);
                    if (j != i && r.get(ix).equals(MI_1)){
                        a = r.plus(c);
                        columns.set(j, a);
                    }
                }
            }
        }
        
        double[][] newMatrix = new double[sm][f];
        for (int i = 0; i < f; i++) {
            DenseVector<ModuloInteger> dv = columns.get(i);
            for (int j = 0; j < sm; j++) {
                newMatrix[j][i] = dv.get(j).equals(MI_1) ? 1 : 0;
            }
        }
        
        for (int d = 0; d < sm; d++) if (!markedRows[d]) dependents.add(d);
        
        return newMatrix;
    }
    
    private BigInteger findFactor(double[][] matrix, int len1, int len2,
        TreeMap<BigInteger, BigInteger> smoothNumbers, BigInteger n) {
        BigInteger factor = BI_ONE, smooth = BI_ONE, product = BI_ONE, sqrtProduct = BI_ONE;
        ArrayList<DenseVector<ModuloInteger>> rows = new ArrayList<DenseVector<ModuloInteger>>();
        ArrayList<Integer> zeroSumRows = new ArrayList();
        int ix = 0, d = 0, numDependents = dependents.size(); boolean factorFound = false;
        ModuloInteger.setModulus(LargeInteger.valueOf(2));
        
        for (int i = 0; i < len1; i++) {
            ArrayList<ModuloInteger> ml = new ArrayList<ModuloInteger>();
            for (int j = 0; j < len2; j++) {
                ml.add(matrix[i][j] == 1 ? MI_1 : MI_0);
            }
            rows.add(DenseVector.valueOf(ml));
        }
        
        while(ix < numDependents && !factorFound){
            d = dependents.get(ix);
            zeroSumRows = getZeroSumRows(d, rows);
            
            smooth = product = sqrtProduct = BI_ONE;
            
            for (int zsr : zeroSumRows) {
                smooth = (BigInteger) smoothNumbers.keySet().toArray()[zsr];
                product = product.multiply(smoothNumbers.get(smooth));
                sqrtProduct = sqrtProduct.multiply(smooth);
            }

            sqrtProduct = sqrtProduct.sqrt();
            factor = gcd(product.subtract(sqrtProduct), n);
            if (!factor.equals(BI_ONE) && !factor.equals(n)) factorFound = true;
            
            ix++;
        }
        
        return factor;
    }
    
    private ArrayList<Integer> getZeroSumRows(int d, ArrayList<DenseVector<ModuloInteger>> rows){
        ArrayList<Integer> zeroSumRows = new ArrayList();
        DenseVector<ModuloInteger> mainRow = rows.get(d), rowToCheck, rowSum, rowTemp;
        
        zeroSumRows.add(d); rowSum = mainRow;
        int numOnes = 0, ix = 0, newNumOnes; 
        for(int i = 0; i < rowSum.getDimension(); i++)
            if(rowSum.get(i).equals(MI_1)) numOnes++;
        
        while(ix < rows.size() && numOnes != 0){
            rowToCheck = rows.get(ix);
            if(ix != d){
                rowTemp = rowToCheck.plus(rowSum);
                newNumOnes = 0;
                for(int i = 0; i < rowTemp.getDimension(); i++)
                    if(rowTemp.get(i).equals(MI_1)) newNumOnes++;
                if (newNumOnes < numOnes) { 
                    zeroSumRows.add(ix);
                    numOnes = newNumOnes;
                    rowSum = rowTemp;
                }
            }
            ix++;
        }
        
        return zeroSumRows;
    }
    
    public record Root (long root, long increment) {
        public Root incrementRoot() { return new Root(root + increment, increment); }
        
        public long getRoot() { return root; }
        
        public long getIncrement() { return increment; }
        
        @Override
        public String toString() { return "[" + root + ", " + increment + "]"; }
    }
}
