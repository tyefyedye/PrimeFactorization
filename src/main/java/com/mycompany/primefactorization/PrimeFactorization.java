package com.mycompany.primefactorization;

/**
 * Cesar Lobaton
 * CSCI 780 - Cryptography
 * RSA Project
 */

import com.mycompany.primefactorization.PollardPMinusOne;
import com.mycompany.primefactorization.QuadraticSieve;
import com.mycompany.primefactorization.PrimeFactoring;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.TreeMap;

public class PrimeFactorization {
    private final static BigInteger BI_ONE = BigInteger.ONE;
    private static ArrayList<BigInteger> numbersToFactor = new ArrayList();
    private static TreeMap<BigInteger, BigInteger> primeFactors = new TreeMap();

    public static void main(String[] args) {
        BigInteger n = new BigInteger(args[0]);
        String method = args[1];
        long bound = Integer.parseInt(args[2]);
        int limit = Integer.parseInt(args[3]);
        BigInteger a, f, f2;

        int count = 0;
        if(method.equals("pollard")) {
            numbersToFactor.add(n);
            PollardPMinusOne pollard = new PollardPMinusOne(bound);
            while((limit != 0 ? count < limit : true) && !numbersToFactor.isEmpty()){
                a = numbersToFactor.removeFirst();
                System.out.printf("Number to factor: %d\n", a);
                f = pollard.getFactor(a);
                f2 = a.divide(f);
                System.out.printf("Factors found: [%d, %d]\n", f, f2);
                checkIfPrimes(f, f2);
                count = (limit != 0 ? count + 1 : 0);
            }
            printPrimeFactors(n);
        } else if (method.equals("quadratic")) {
            numbersToFactor.add(n);
            QuadraticSieve q = new QuadraticSieve();
            while((limit != 0 ? count < limit : true) && !numbersToFactor.isEmpty()){
                a = numbersToFactor.removeFirst();
                System.out.printf("Number to factor: %d\n", a);
                f = q.getFactor(a, bound);
                f2 = a.divide(f);
                System.out.printf("Factors found: [%d, %d]\n", f, f2);
                checkIfPrimes(f, f2);
                count = (limit != 0 ? count + 1 : 0);
            }
            printPrimeFactors(n);
        } 
    }
    
    public static void checkIfPrimes(BigInteger f, BigInteger f2){
        PrimeFactoring pf = new PrimeFactoring();
        boolean isPrime = false, isPrime2 = false;
        isPrime = pf.checkIfPrime(f, 100000);
        isPrime2 = pf.checkIfPrime(f2, 100000);
        if (isPrime){
            System.out.printf("%d is prime\n", f);
            if (primeFactors.containsKey(f)) 
                primeFactors.put(f, primeFactors.get(f).add(BI_ONE));
            else
                primeFactors.put(f, BI_ONE);
        } else {
            numbersToFactor.add(f);
        }
        if (isPrime2){
            System.out.printf("%d is prime\n", f2);
            if (primeFactors.containsKey(f2)) 
                primeFactors.put(f, primeFactors.get(f2).add(BI_ONE));
            else
                primeFactors.put(f2, BI_ONE);
        } else {
            numbersToFactor.add(f2);
        }
    }
    
    public static void printPrimeFactors(BigInteger n){
        BigInteger fPow;
        System.out.printf("Factors of %d are: ", n);
        for (BigInteger factor : primeFactors.keySet()){
            fPow = primeFactors.get(factor);
            System.out.print(factor + (fPow.equals(BI_ONE) ? "" : "^" + fPow) + 
            (factor.equals(primeFactors.lastKey()) ? "" : " * "));
        }
    }
}
