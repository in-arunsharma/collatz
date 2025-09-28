#include <iostream>
#include <set>
using namespace std;

int collatz_steps(long long original_n) {
    long long n = original_n;
    int steps = 0;
    set<long long> seen_numbers;  // Track numbers we've seen
    
    while (n != 1) {
        // Check if we've seen this number before (cycle detection!)
        if (seen_numbers.find(n) != seen_numbers.end()) {
            cout << "CYCLE DETECTED! Starting number: " << original_n << endl;
            cout << "Cycle involves number: " << n << " after " << steps << " steps" << endl;
            cout << "This could be a counterexample to Collatz Conjecture!" << endl;
            return -1;  // Special return value for cycles
        }
        
        seen_numbers.insert(n);  // Remember this number
        
        if (n % 2 == 0) {
            n = n / 2;  // Even: divide by 2
        } else {
            n = 3 * n + 1;  // Odd: multiply by 3 and add 1
        }
        steps++;
        
    }
    
    return steps;
}

int main() {
    cout << "Basic Collatz Conjecture Calculator with Cycle Detection" << endl;
    cout << "=======================================================" << endl;
    
    long long counterexamples_found = 0;
    long long numbers_tested = 0;
    
    // Test with some small numbers first  
    for (long long i = 1000000; i <= 2000000; i++) {  // Reduced for testing with cycle detection
        int steps = collatz_steps(i);
        numbers_tested++;
        
        if (steps == -1) {  // Cycle detected!
            counterexamples_found++;
            cout << "*** POTENTIAL COUNTEREXAMPLE FOUND ***" << endl;
        } 

        
        // Progress indicator
        //if (i % 10000 == 0) {
        //    cout << "Tested " << i << " numbers so far..." << endl;
        //}
    }
    
    cout << endl;
    cout << "=== FINAL RESULTS ===" << endl;
    cout << "Numbers tested: " << numbers_tested << endl;
    cout << "Counterexamples found: " << counterexamples_found << endl;
    
    if (counterexamples_found == 0) {
        cout << "No counterexamples found - Collatz Conjecture holds for tested range!" << endl;
    } else {
        cout << "AMAZING! Found " << counterexamples_found << " counterexamples!" << endl;
    }
       
    return 0;
}
