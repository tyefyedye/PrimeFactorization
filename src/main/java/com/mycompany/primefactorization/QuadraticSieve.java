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
    
    private ArrayList<Integer> dependents = new ArrayList();
    private long lastBoundUsed;
    
    public QuadraticSieve(){ }
    
    public BigInteger getFactor(BigInteger n, long b, BigInteger lsb, 
            BigInteger usb, boolean useShanksTonelli, int interval){
        if(n.equals(BI_NINE)) return BI_THREE;
        else if (n.compareTo(BI_TEN) == -1) return BI_TWO;
        
        long bound; ArrayList<Long> primes;
        BigInteger p, pMinusOne, qr, a, factor;
        ArrayList<BigInteger> factorBase = new ArrayList();
        ArrayList<Root> roots = new ArrayList();
        TreeMap<BigInteger, BigInteger> smoothNumbers = new TreeMap();
        TreeMap<BigInteger, double[]> smoothNumbersPrimeFactors = new TreeMap();
        
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
                if (qr.equals(BI_ONE) || qr.equals(BI_ZERO)) factorBase.add(p);
            }
            numFactors = factorBase.size();
            System.out.println(numFactors + " quadratic residues found.");
            if (numFactors < 7) {
                System.out.println("Not enough quadratic residues found.");
                factorBase.clear(); bound = bound + 10;
                System.out.println("Trying bound B = " + bound);
            }
        }
//        for(BigInteger f : factorBase) System.out.print(f + " ");
//        System.out.println();
        
        a = n.sqrt();
//        System.out.println(a);

        if (useShanksTonelli) {
            System.out.println("Getting roots over factor base...");
            BigInteger[] squareRoots;
            BigInteger r1, r2;
            for (BigInteger fb : factorBase){
                if (fb.equals(BI_TWO)) roots.add(new Root(1, 2));
                else if(n.mod(fb).equals(BI_ZERO)) 
                    roots.add(new Root(0,fb.longValue()));
                else {
                    squareRoots = tonelliShanks(n, fb);
                    r1 = squareRoots[0].subtract(a).mod(fb).add(lsb);
                    r2 = squareRoots[1].subtract(a).mod(fb).add(lsb);
                    roots.add(new Root(r1.longValue(), fb.longValue()));
                    roots.add(new Root(r2.longValue(), fb.longValue()));
                }
            }
            System.out.println("Roots found");
        }
        
        System.out.println("Finding smooth numbers...");
        int numSmoothNumbers = 0, i = 0, j = 0, c = 0;
        boolean allRootsExceedUpperBound = false;
        BigInteger root, f, t, r, sb = lsb;
        double[] pExp; 
        while ((useShanksTonelli ? !allRootsExceedUpperBound : !sb.equals(usb)) && numSmoothNumbers < numFactors){
            root = (useShanksTonelli ? BigInteger.valueOf(roots.get(i % numFactors).getRoot()) : sb);
            t = polynomial(a, root, n);
            r = t;
            if(root.compareTo(usb) != 1 && !smoothNumbers.containsKey(t)){
                c++;
                j = 0;
                pExp = new double[numFactors];
                Arrays.fill(pExp, 0);
                while (j < numFactors && !r.equals(BI_ZERO) && !r.equals(BI_ONE)){
                    f = factorBase.get(j);
                    if (!r.mod(f).equals(BI_ZERO)) j++;
                    else { r = r.divide(f); pExp[j] = (pExp[j] + 1) % 2; }
                }
                if (j < numFactors){
                    smoothNumbers.put(t, a.add(root));
                    smoothNumbersPrimeFactors.put(t, pExp);
                    numSmoothNumbers++;
                    if (numSmoothNumbers % interval == 0)
                        System.out.println(numSmoothNumbers + " smooth numbers found");
                }
            }
            if (useShanksTonelli){
                roots.set(i % numFactors, roots.get(i % numFactors).incrementRoot());
                i++;
                if (c == 0) allRootsExceedUpperBound = true; 
                else if (i % numFactors == 0) c = 0;
            } else sb = sb.add(BI_ONE);
        }
        System.out.println(numSmoothNumbers + " total smooth numbers found");
        
        System.out.println("Building matrix...");
        int k = 0; double[] arr;
        double[][] matrix = new double[numFactors][numSmoothNumbers];
        for(BigInteger sn : smoothNumbers.keySet()){
            arr = smoothNumbersPrimeFactors.get(sn);
            if (arr != null) { /*System.out.println(sn + " " + smoothNumbers.get(sn));*/ matrix[k] = arr; k++; } 
        }
        System.out.println("Matrix built");
        
//        for (int g = 0; g < numFactors; g++){
//            for (int h = 0; h < numFactors; h++){
//                System.out.print((int) matrix[g][h] + " ");
//            }
//            System.out.println();
//        }
        
        System.out.println("Performing Gaussian Elimination on matrix...");
        dependents.clear();
        matrix = gaussianElimination(matrix, numFactors, numSmoothNumbers);
        System.out.println("Gaussian Elimination complete");
        
//        for (int g = 0; g < numFactors; g++){
//            for (int h = 0; h < numFactors; h++){
//                System.out.print((int) matrix[g][h] + " ");
//            }
//            System.out.println();
//        }
//        
//        for(int d : dependents) System.out.print(d + " ");
//        System.out.println();
        
        System.out.println("Finding factor...");
        factor = findFactor(matrix, numFactors, numSmoothNumbers, smoothNumbers, n);
        if(!factor.equals(BI_ONE) && !factor.equals(n)){
            System.out.println(factor);
            System.out.println("Factor found!");
        }
        
        return factor;
    }
    
    public long incrementBound() { return lastBoundUsed + 10; }
    
    private BigInteger polynomial(BigInteger a, BigInteger i, BigInteger n) {
        BigInteger ai = a.add(i);
        return ai.pow(2).subtract(n);
    }
    
    private double[][] gaussianElimination(double[][] matrix, int len1, int len2){
        ModuloInteger.setModulus(LargeInteger.valueOf(2));
        int ix; DenseVector<ModuloInteger> c, r, a;
        ArrayList<DenseVector<ModuloInteger>> columns = new ArrayList();
        boolean markedRows[] = new boolean[len1]; Arrays.fill(markedRows, false);

        for (int i = 0; i < len2; i++) {
            ArrayList<ModuloInteger> ml = new ArrayList<ModuloInteger>();
            for (int j = 0; j < len1; j++) {
                ml.add(matrix[j][i] == 1 ? MI_1 : MI_0);
            }
            columns.add(DenseVector.valueOf(ml));
        }
        
        for (int i = 0; i < len2; i++){
            ix = 0; c = columns.get(i);
            while (ix < len2 && !c.get(ix).equals(MI_1)) ix++;
            if (ix < len2) {
                markedRows[ix] = true;
                for (int j = 0; j < columns.size(); j++){
                    r = columns.get(j);
                    if (!r.equals(c) && r.get(ix).equals(MI_1)){
                        a = r.plus(c);
                        columns.set(j, a);
                    }
                }
            }
        }
        
        double[][] newMatrix = new double[len1][len2];
        for (int i = 0; i < len2; i++) {
            DenseVector<ModuloInteger> dv = columns.get(i);
            for (int j = 0; j < len1; j++) {
                newMatrix[j][i] = dv.get(j).equals(MI_1) ? 1 : 0;
            }
            if (!markedRows[i]) dependents.add(i);
        }
        
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
            
//            for(int zsr : zeroSumRows) System.out.print(zsr + " ");
//            System.out.println();
            
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
