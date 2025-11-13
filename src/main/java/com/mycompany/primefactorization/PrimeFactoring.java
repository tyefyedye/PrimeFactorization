package com.mycompany.primefactorization;

import java.math.BigInteger;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Random;

public abstract class PrimeFactoring {
    public final BigInteger BI_ZERO = BigInteger.ZERO;
    public final BigInteger BI_ONE = BigInteger.ONE;
    public final BigInteger BI_TWO = BigInteger.TWO;
    public final BigDecimal BI_TEN = BigDecimal.valueOf(10);
    
    public BigInteger gcd(BigInteger n, BigInteger m){
        if (n.equals(BI_ZERO)) return m;
        BigInteger a = n.max(m), b = n.min(m), r = a.mod(b);
        
        while (r.compareTo(BI_ZERO) == 1) {
            a = b; b = r; r = a.mod(b);
        }
        return b;
    }
    
    public ArrayList<Long> sieveOfEratosthenes(long b){
        ArrayList<Long> primes = new ArrayList();
        boolean[] numsToCheck = new boolean[(int) b];
        Arrays.fill(numsToCheck, true);
        
        //for(long i = 2; i <= b; i++) primes.add(i);

        int limit = (int) Math.round(Math.sqrt(b));
        
        int p;
        for (int i = 2; i <= limit; i++){
            if(numsToCheck[i]) {
                p = (int) Math.pow(i, 2);
                for (int j = p; j < b; j++)
                    if (j % i == 0) numsToCheck[j] = false;
            }
        }
        
        for (int k = 2; k < b; k++) if(numsToCheck[k]) primes.add((long) k);
        
//        int i = 0; long p, q;
//        while(primes.get(i) <= limit){
//            p = primes.get(i);
//            for(int j = 0; j < primes.size(); j++){
//                q = primes.get(j);
//                if (q != p && q % p == 0) primes.remove(j);
//            }
//            i++;
//        }
        return primes;
    }
    
    public BigInteger[] tonelliShanks(BigInteger n, BigInteger p){
        BigInteger[] squareRoots = new BigInteger[2];
        BigInteger pMinusOne = p.subtract(BI_ONE), t = BI_ZERO, y = BI_TWO,
            x = BI_ZERO, z = BI_ZERO, b = BI_ZERO, w = BI_ZERO;
        int s = 1, e, q;
        boolean tIsOdd = false;
        
        while(BI_TWO.pow(s).compareTo(pMinusOne) != 1 && !tIsOdd) {
            t = pMinusOne.divide(BI_TWO.pow(s));
            if (t.mod(BI_TWO).equals(BI_ONE)) tIsOdd = true;
            else s++;
        }

        while (y.modPow(pMinusOne.divide(BI_TWO), p).equals(BI_ONE))
            y = y.add(BI_ONE);
        
        x = n.modPow(t.add(BI_ONE).divide(BI_TWO), p);
        z = y.modPow(t, p);
        b = n.modPow(t, p);
        e = s;
        
        while(!b.equals(BI_ONE)){
            q = 1;
            while (!b.modPow(BI_TWO.pow(q), p).equals(BI_ONE)) q++;
            w = z.modPow(BI_TWO.pow(e-q-1), p);
            x = x.multiply(w).mod(p);
            b = b.multiply(w.modPow(BI_TWO, p)).mod(p);
            z = w.modPow(BI_TWO, p);
            e = q;
        }
        
        squareRoots[0] = x;
        squareRoots[1] = p.subtract(x);
        
        return squareRoots;
    }
    
    public BigInteger getRandomNum(BigInteger n){
        int bitLength = n.bitLength();
        Random rand = new Random();
        BigInteger randNum = new BigInteger(bitLength, rand);
        return randNum;
    }
    
    public long getBound(BigInteger n){
        int numDigits = n.toString().length();
        
        BigDecimal bd = new BigDecimal(n);
        double val = bd.divide(BI_TEN.pow(numDigits)).doubleValue();
        double logN = numDigits * Math.log(10) + Math.log(val);
        
        double l = Math.exp(Math.sqrt(logN * Math.log(logN)));
        
        return (long) Math.ceil(Math.pow(l, 1/Math.sqrt(2)));
    }
    
    public boolean checkIfPrime(BigInteger n, int bound) {
        ArrayList<Long> primes = sieveOfEratosthenes(bound);
        int numPrimes = primes.size(), s = 1;
        boolean mIsOdd = false, testResult = false;
        BigInteger maxPrime = BigInteger.valueOf(primes.get(numPrimes-1));
        BigInteger p, nMinusOne, m = BI_ZERO;
        
        if (n.compareTo(maxPrime) == -1)
            for (long prime : primes){
                p = BigInteger.valueOf(prime);
                if (n.equals(p)) return true;
            }
        
        nMinusOne = n.subtract(BI_ONE);
        while(BI_TWO.pow(s).compareTo(nMinusOne) != 1 && !mIsOdd) {
            m = nMinusOne.divide(BI_TWO.pow(s));
            if (m.mod(BI_TWO).equals(BI_ONE)) mIsOdd = true;
            else s++;
        }
        
        for (long prime : primes){
            p = BigInteger.valueOf(prime);
            testResult = rabinMillerPrimalityTest(n, p, nMinusOne, m, s);
            if (!testResult) return false;
        }
        return true;
    }
    
    public boolean rabinMillerPrimalityTest(BigInteger n, BigInteger a,
        BigInteger nMinusOne, BigInteger m, int s){
        BigInteger c;
        boolean isProbablePrime = false;

        c = a.modPow(m, n);
        if (c.equals(BI_ONE) || c.equals(nMinusOne)) isProbablePrime = true;
        else {
            for (int i = 0; i < s; i++){
                c = c.pow(2).mod(n);
                if (c.equals(nMinusOne)) { isProbablePrime = true; break; }
            }
        }
        return isProbablePrime;
    }
}
