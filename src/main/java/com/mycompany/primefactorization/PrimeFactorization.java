package com.mycompany.primefactorization;

/**
 * Cesar Lobaton
 * CSCI 780 - Cryptography
 * RSA Project
 */

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.Scanner;

public class PrimeFactorization {
    private final static BigInteger BI_ONE = BigInteger.ONE;
    private static ArrayList<BigInteger> numbersToFactor = new ArrayList();
    private static TreeMap<BigInteger, BigInteger> primeFactors = new TreeMap();

    public static void main(String[] args) {
        BigInteger n = new BigInteger(args[0]);
        String method = args[1];
        long bound = Integer.parseInt(args[2]);
        int limit = Integer.parseInt(args[3]);
        int range = Integer.parseInt(args[4]);
        int threshold = Integer.parseInt(args[5]);
        int minCandidates = Integer.parseInt(args[6]);
        int reductionRate = Integer.parseInt(args[7]);
        BigInteger a, f, f2;

        int count = 0; String whatToIncrement;
        Scanner input = new Scanner(System.in);
        if(method.equals("pollard")) {
            numbersToFactor.add(n);
            PollardPMinusOne pollard = new PollardPMinusOne(bound);
            while((limit != 0 ? count < limit : true) && !numbersToFactor.isEmpty()){
                a = numbersToFactor.removeFirst();
                System.out.printf("Number to factor: %d\n", a);
                f = pollard.getFactor(a);
                f2 = a.divide(f);
                System.out.printf("Returned factors: [%d, %d]\n", f, f2);
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
                f = q.getFactor(a, bound, range, threshold, minCandidates, reductionRate);
                f2 = a.divide(f);
                if (f.equals(BI_ONE) || f2.equals(BI_ONE)){
                    System.out.print("No factors found. Increment bound range, threshold, or candidates?? ");
                    whatToIncrement = input.nextLine();
                    if (whatToIncrement.equals("bound")){
                        System.out.print("Input new bound: ");
                        bound = input.nextLong();
                    } else if (whatToIncrement.equals("range")){
                        System.out.print("Input new range: ");
                        range = input.nextInt();
                    } else if (whatToIncrement.equals("threshold")){
                        System.out.print("Input new threshold: ");
                        threshold = input.nextInt();
                    } else if (whatToIncrement.equals("candidates")){
                        System.out.print("Input new candidates: ");
                        minCandidates = input.nextInt();
                    }
                    numbersToFactor.add(a);
                } else {
                    bound = Integer.parseInt(args[2]);
                    System.out.printf("Factors found: [%d, %d]\n", f, f2);
                    checkIfPrimes(f, f2);
                    count = (limit != 0 ? count + 1 : 0);
                }
            }
            input.close();
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
                primeFactors.put(f2, primeFactors.get(f2).add(BI_ONE));
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
