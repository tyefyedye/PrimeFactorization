
package com.mycompany.primefactorization;

/**
 * Cesar Lobaton
 * CSCI 780 - Cryptography
 * RSA Project
 */

import java.math.BigInteger;
import java.util.ArrayList;

public class PollardPMinusOne extends PrimeFactoring {
    private long bound;
    
    public PollardPMinusOne(long bound){
        this.bound = bound;
    }
    
    public BigInteger getFactor(BigInteger n){
        BigInteger a = BI_TWO, aM = BI_ZERO, d = BI_ZERO, e = BI_ONE, m, p;
        ArrayList<Long> primes;
        long b = bound; int exp;
        
        while(e.equals(BI_ONE) || e.equals(n)) {
            System.out.printf("Try a = %d, b = %d\n", a, b);
            primes = sieveOfEratosthenes(b);
            m = BI_ONE;
            for (int i = 0; i < primes.size(); i++){
                exp = (int) Math.floor(Math.log(b) / Math.log(primes.get(i)));
                p = BigInteger.valueOf(primes.get(i));
                m = m.multiply(p.pow(exp));
            }

            d = gcd(a,n);
            if (!d.equals(BI_ONE)) return d;
            else {
                aM = a.modPow(m, n);
                e = gcd(aM.subtract(BI_ONE), n);
                if(e.equals(BI_ONE)) b += b;
                else if(e.equals(n)) { a = getRandomNum(n); b = bound; }
            }
        }
        return e;
    }
}
