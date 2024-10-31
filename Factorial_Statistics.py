"""
Factorial Product Search Program - Revision 1.4
=============================================
This program searches for solutions to the equation a! x b! = c! where 1 < a < b < c.
Uses prime factor exponent lists and sliding window technique to manage memory usage.
Collects statistics to demonstrate vanishing possibility of solutions.

Example: 6! = 720 = 2^4 x 3^2 x 5^1 is represented as [4, 2, 1]
- 4 is the exponent of the first prime (2)
- 2 is the exponent of the second prime (3)
- 1 is the exponent of the third prime (5)

The only known non-trivial solution is 6! x 7! = 10!

Author: Ken Clements
Date: October 30, 2024
"""

import math
import csv
from datetime import datetime
import argparse

# Search Parameters
SEARCH_LIMIT = 400  # Upper limit for factorial search
WINDOW_SIZE = 10000  # Size of sliding window for factorial lists
WINDOW_OVERLAP = 1000  # Number of entries to overlap between windows
MAX_B_LOOKBACK = 1000  # How far back to look for b values from c
MAX_A_LOOKBACK = 100  # How far back to look for a values from b

# Progress Reporting
FACTOR_PROGRESS_INTERVAL = 1000  # How often to report factorial generation progress

# Statistics Collection
STATS_ENABLED = True  # Whether to collect and save statistics
STATS_FILE_PREFIX = "factorial_search_stats"  # Prefix for stats file names

class SearchStats:
    def __init__(self):
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.filename = f"{STATS_FILE_PREFIX}_{timestamp}.csv"
        self.current_file = None
        self.headers = ['c', 'max_prime_idx', 'b_range_size',
                       'length_mismatches', 'internal_zeros',
                       'last_exp_too_large', 'b_attempts',
                       'smaller_rejects', 'avg_ba_diff']
        
        # Create file and write headers
        with open(self.filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(self.headers)

    def open_stats_file(self):
        """Open the stats file in append mode"""
        if self.current_file is None:
            self.current_file = open(self.filename, 'a', newline='')
            self.writer = csv.writer(self.current_file)

    def close_stats_file(self):
        """Close the stats file if it's open"""
        if self.current_file is not None:
            self.current_file.close()
            self.current_file = None
            self.writer = None

    def record_c_stats(self, c, stats):
        """Record statistics for a single c value"""
        if self.current_file is None:
            self.open_stats_file()
        avg_ba_diff = (sum(stats['ba_diffs']) / len(stats['ba_diffs'])
                      if stats['ba_diffs'] else 0)
        self.writer.writerow([
            c, stats['max_prime_idx'], stats['b_range_size'],
            stats['length_mismatches'], stats['internal_zeros'],
            stats['last_exp_too_large'], stats['b_attempts'],
            stats['smaller_rejects'], avg_ba_diff
        ])

    def __del__(self):
        """Ensure file is closed when object is destroyed"""
        self.close_stats_file()

class FactorialSearch:
    def __init__(self, limit: int = SEARCH_LIMIT):
        self.limit = limit
        self.factorial_primes_list = [[0], [0], [1]]  # Start with 2! at index 2
        self.prime_cache = [2]  # Initialize with first prime
        self.window_start = 0  # Offset for indexing
        self.stats = SearchStats() if STATS_ENABLED else None

    def get_nth_prime(self, nth: int) -> int:
        """Get the nth prime number (0-based indexing)."""
        if nth < len(self.prime_cache):
            return self.prime_cache[nth]
        
        size_factor = 2
        s = (nth * size_factor)
        
        def get_primes(limit):
            sieve = [True] * limit
            for i in range(2, int(limit ** 0.5) + 1):
                if sieve[i]:
                    for j in range(i*i, limit, i):
                        sieve[j] = False
            return sieve
        
        while True:
            primes = get_primes(s)
            new_primes = [i for i in range(2, s) if primes[i]]
            if len(new_primes) > nth:
                self.prime_cache.extend(new_primes[len(self.prime_cache):])
                return self.prime_cache[nth]
            size_factor += 1
            s = (nth * size_factor)

    def get_p_factor_list(self, n: int) -> list:
        """Get the prime factor exponent list for n."""
        if n < 2:
            return []
        
        result = []
        prime_index = 0
        working_n = n
        max_prime = 0
        
        # First find the largest prime factor
        temp_n = n
        while temp_n > 1:
            prime = self.get_nth_prime(prime_index)
            if temp_n % prime == 0:
                max_prime = prime_index
                while temp_n % prime == 0:
                    temp_n //= prime
            prime_index += 1
        
        # Build prime factor list
        for i in range(max_prime + 1):
            prime = self.get_nth_prime(i)
            count = 0
            while working_n % prime == 0:
                count += 1
                working_n //= prime
            result.append(count)
        
        return result

    def get_factorial_factors(self, idx):
        """Get factorial factors adjusting for window offset"""
        return self.factorial_primes_list[idx - self.window_start]

    def slide_window(self, new_start, new_end):
        """Slide the window of factorial prime factors"""
        keep_start = max(2, new_start - WINDOW_OVERLAP)
        if keep_start > self.window_start:
            remove_count = keep_start - self.window_start
            self.factorial_primes_list = self.factorial_primes_list[remove_count:]
            self.window_start = keep_start

        current_end = self.window_start + len(self.factorial_primes_list)
        
        for k in range(current_end, new_end + 1):
            k_factors = self.get_p_factor_list(k)
            prev_factorial_factors = self.get_factorial_factors(k-1).copy()
            new_factorial_factors = sum_factors(k_factors, prev_factorial_factors)
            self.factorial_primes_list.append(new_factorial_factors)
            
            if k % FACTOR_PROGRESS_INTERVAL == 0:
                print(f"Generated factorial prime factors up to {k}")

    def search(self):
        """Main search function with sliding window"""
        window_size = WINDOW_SIZE
        current_start = 2
        total_exceptions = 0
        exceptions_past_10 = 0
        
        while current_start < self.limit:
            current_end = min(current_start + window_size, self.limit)
            print(f"\nSliding window to range {current_start}! to {current_end}!")
            
            self.slide_window(current_start, current_end)
            
            if self.stats:
                self.stats.open_stats_file()
            
            exceptions = search_section(self, current_start, current_end)
            
            if exceptions:
                for a, b, c in exceptions:
                    print(f"Found exception: {a}! × {b}! = {c}!")
                    total_exceptions += 1
                    if c > 10:
                        exceptions_past_10 += 1
            else:
                print("No exceptions found in this range")
            
            if self.stats:
                self.stats.close_stats_file()
            
            if current_end >= self.limit:
                break
            
            current_start = current_end - WINDOW_OVERLAP
        
        print("\nSearch Complete!")
        print(f"Searched factorial products up to {self.limit}!")
        print(f"Total exceptions found: {total_exceptions}")
        print(f"Exceptions found past 10!: {exceptions_past_10}")

def remove_trailing_zeros(diff):
    while diff and diff[-1] == 0:
        diff.pop()
    return diff

def sum_factors(list1, list2):
    list1 = list1.copy()
    list2 = list2.copy()
    length = max(len(list1), len(list2))
    list1 += [0] * (length - len(list1))
    list2 += [0] * (length - len(list2))
    result = [a + b for a, b in zip(list1, list2)]
    while result and result[-1] == 0:
        result.pop()
    return result

def search_section(searcher, start_idx, end_idx):
    exceptions = []
    for c_idx in range(start_idx, end_idx):
        # Initialize statistics for this c value
        stats_for_c = {
            'max_prime_idx': 0,  # Highest prime index in c's factors
            'b_range_size': 0,  # Size of b range for this c
            'length_mismatches': 0,  # Count of b values rejected due to length
            'internal_zeros': 0,  # Count of diff lists with internal zeros
            'last_exp_too_large': 0,  # Count of diff lists with last exp > 1
            'b_attempts': 0,  # Count of b values that passed length check
            'smaller_rejects': 0,  # Count of a values rejected for small exponents
            'ba_diffs': []  # Store b-a differences when rejecting
        }
        
        c_factors = searcher.get_factorial_factors(c_idx)
        
        # Find highest non-zero prime index in c_factors
        for i in range(len(c_factors)-1, -1, -1):
            if c_factors[i] != 0:
                stats_for_c['max_prime_idx'] = i
                break
        
        # Calculate b range size
        first_b = c_idx - 2
        last_b = max(2, c_idx-MAX_B_LOOKBACK)
        stats_for_c['b_range_size'] = first_b - last_b + 1
        
        for b_idx in range(first_b, last_b, -1):
            b_factors = searcher.get_factorial_factors(b_idx)
            
            if len(b_factors) != len(c_factors):
                stats_for_c['length_mismatches'] += 1
                continue
            
            stats_for_c['b_attempts'] += 1
            diff = [c_factors[j] - b_factors[j] for j in range(len(c_factors))]
            diff = remove_trailing_zeros(diff.copy())
            
            if any(x == 0 for x in diff[:-1]):
                stats_for_c['internal_zeros'] += 1
                continue
            
            if diff and diff[-1] > 1:
                stats_for_c['last_exp_too_large'] += 1
                continue
            
            # Search backwards from b-1 for matching a
            found_smaller = False
            for a_idx in range(b_idx-1, max(2, b_idx-MAX_A_LOOKBACK), -1):
                a_factors = searcher.get_factorial_factors(a_idx)
                
                if len(a_factors) < len(diff):
                    found_smaller = True
                    stats_for_c['smaller_rejects'] += 1
                    stats_for_c['ba_diffs'].append(b_idx - a_idx)
                    break
                
                if len(a_factors) > len(diff):
                    continue
                
                # Check each exponent in order
                for i in range(len(diff)):
                    if a_factors[i] < diff[i]:
                        found_smaller = True
                        stats_for_c['smaller_rejects'] += 1
                        stats_for_c['ba_diffs'].append(b_idx - a_idx)
                        break
                    if a_factors[i] > diff[i]:
                        break
                else:  # All exponents match
                    exceptions.append((a_idx, b_idx, c_idx))
                    break
                
                if found_smaller:
                    break
        
        # Record stats for this c value if we had any activity
        if searcher.stats and (stats_for_c['b_attempts'] > 0 or
                             stats_for_c['length_mismatches'] > 0):
            searcher.stats.record_c_stats(c_idx, stats_for_c)
    
    return exceptions

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Factorial Product Search Program - Searches for solutions to a! × b! = c!',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--limit', type=int, default=SEARCH_LIMIT,
                       help='Upper limit for factorial search')
    parser.add_argument('--window-size', type=int, default=WINDOW_SIZE,
                       help='Size of sliding window for factorial lists')
    parser.add_argument('--no-stats', action='store_true',
                       help='Disable statistics collection')
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    STATS_ENABLED = not args.no_stats
    searcher = FactorialSearch(limit=args.limit)
    searcher.search()

# END OF PROGRAM

