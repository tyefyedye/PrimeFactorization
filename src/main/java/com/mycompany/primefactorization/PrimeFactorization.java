package com.mycompany.primefactorization;

/**
 * Cesar Lobaton
 * CSCI 780 - Cryptography
 * RSA Project
 */

import com.mycompany.primefactorization.PollardPMinusOne;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.TreeMap;

public class PrimeFactorization {
    private final static BigInteger BI_ONE = BigInteger.ONE;

    public static void main(String[] args) {
        BigInteger n = new BigInteger(args[0]);
        String method = args[1];
        long bound = Integer.parseInt(args[2]);

        if(method.equals("pollard")) {
            BigInteger a, f, f2; boolean isPrime, isPrime2;
            ArrayList<BigInteger> numbersToFactor = new ArrayList();
            TreeMap<BigInteger, BigInteger> primeFactors = new TreeMap();
            
            numbersToFactor.add(n);
            while(!numbersToFactor.isEmpty()){
                PollardPMinusOne pollard = new PollardPMinusOne(bound);
                a = numbersToFactor.removeFirst();
                System.out.printf("Number to factor: %d\n", a);
                f = pollard.getFactor(a);
                f2 = a.divide(f);
                System.out.printf("Factors found: [%d, %d]\n", f, f2);
                
                isPrime = isPrime2 = false;
                isPrime = pollard.checkIfPrime(f, 1000);
                isPrime2 = pollard.checkIfPrime(f2, 1000);
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
            BigInteger fPow;
            System.out.printf("Factors of %d are: ", n);
            for (BigInteger factor : primeFactors.keySet()){
                fPow = primeFactors.get(factor);
                System.out.print(factor + (fPow.equals(BI_ONE) ? "" : "^" + fPow) + 
                (factor.equals(primeFactors.lastKey()) ? "" : " * "));
            }
        } else if (method.equals("quadratic")) {
            System.out.printf("Number to factor: %d\n", n);
            QuadraticSieve q = new QuadraticSieve();
            BigInteger a = q.getFactor(bound, n);
        } 
    }
}
