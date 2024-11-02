"""
Factorial Product Search Program - Revision 3.1
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
Date: Nov. 1, 2024
"""

from tqdm import tqdm
import math
import csv
from datetime import datetime
import argparse

# Search Parameters
SEARCH_LIMIT = 2500  # Upper limit for factorial search
WINDOW_SIZE = 50000  # Size of sliding window for factorial lists
WINDOW_OVERLAP = 5000  # Number of entries to overlap between windows
MAX_B_LOOKBACK = 1000  # How far back to look for b values from c


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
        self.headers = ['c', 'c_list_length',
                       'length_mismatches', 'non-monotonic', 'internal_zeros',
                       'last_exp_too_large', 'b_attempts',
                       'a_attempts',
                       'smaller_rejects', 'PF-slope']
        
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
        
        self.writer.writerow([
            c, stats['c_list_length'],
            stats['length_mismatches'], stats['non-monotonic'], stats['internal_zeros'],
            stats['last_exp_too_large'], stats['b_attempts'],
            stats['a_attempts'],
            stats['smaller_rejects'], stats['PF-slope']
        ])

    def __del__(self):
        """Ensure file is closed when object is destroyed"""
        self.close_stats_file()

# END OF CLASS SearchStats

class FactorialSearch:
    def __init__(self, limit: int = SEARCH_LIMIT, use_progress_bar: bool = False, k_primes: int = 0):
        self.limit = limit
        self.factorial_primes_list = [[0], [0], [1]]  # Start with 2! at index 2
        self.prime_cache = [2]  # Initialize with first prime
        self.window_start = 0  # Offset for indexing
        self.stats = SearchStats() if STATS_ENABLED else None
        self.use_progress_bar = use_progress_bar
        self.k_primes = k_primes


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


    def slide_window(self, new_start, new_end, quiet=False):  # Add quiet parameter
        """Slide the window of factorial prime factors"""
        keep_start = max(2, new_start - WINDOW_OVERLAP)
        if keep_start > self.window_start:
            remove_count = keep_start - self.window_start
            self.factorial_primes_list = self.factorial_primes_list[remove_count:]
            self.window_start = keep_start

        current_end = self.window_start + len(self.factorial_primes_list)
    
        # Set up progress tracking for factorial generation
        remaining_factors = range(current_end, new_end + 1)
        if self.use_progress_bar:
           progress_iter = tqdm(remaining_factors, 
                           desc="Generating factorial prime factors",
                           unit="factorial")
        else:
            progress_iter = remaining_factors
    
        for k in progress_iter:
            k_factors = self.get_p_factor_list(k)
            prev_factorial_factors = self.get_factorial_factors(k-1).copy()
            new_factorial_factors = sum_factors(k_factors, prev_factorial_factors)
            self.factorial_primes_list.append(new_factorial_factors)
        
            if not self.use_progress_bar and k % FACTOR_PROGRESS_INTERVAL == 0 and not quiet:
                # Find the index of the last non-zero exponent
                last_prime_idx = len(new_factorial_factors) - 1
                highest_prime = self.prime_cache[last_prime_idx]
                print(f"Found prime factors up to {k}! Highest prime used is {highest_prime}")

    def search(self, quiet=False):  # Add quiet parameter
        """Main search function with sliding window"""
        window_size = WINDOW_SIZE
        current_start = 2
        total_exceptions = 0
        exceptions_past_10 = 0
    
        while current_start < self.limit:
            current_end = min(current_start + window_size, self.limit)
            if not quiet:
                print(f"\nProcessing window {current_start-1}! to {current_end}!")

        
            self.slide_window(current_start, current_end)
        
            if self.stats:
                self.stats.open_stats_file()
        
            exceptions = search_section(self, current_start, current_end)
        
            if exceptions:
                for a, b, c in exceptions:
                    if not quiet:
                        if self.k_primes == 0:
                            print(f"\nFound exception: {a}! x {b}! = {c}!")
                        else:
                            print(f"Found possible exception: {a}! x {b}! = {c}!")

                    total_exceptions += 1
                    if c > 10:
                        exceptions_past_10 += 1
        
            if self.stats:
                self.stats.close_stats_file()
        
            if current_end >= self.limit:
                break
        
            current_start = current_end - WINDOW_OVERLAP
        if not quiet:
            print("\nSearch Complete!")
            print(f"Searched factorial products up to {self.limit-1}!")
            print(f"Total exceptions found: {total_exceptions}")
            print(f"Exceptions found past 10!: {exceptions_past_10}")


# END OF CLASS FactorialSearch


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

def search_section(searcher, start_idx, end_idx, quiet=False):
    exceptions = []
    
    # Set up progress tracking for search phase
    if searcher.use_progress_bar:
        c_range = tqdm(range(start_idx, end_idx),
                      desc="Searching for solutions",
                      unit="c")
    else:
        c_range = range(start_idx, end_idx)
     
    k_primes = searcher.k_primes # if > 0, only searh for match in lower prime factors
    # print(f"k_primes = {k_primes}")

    for c_idx in c_range:
        # Initialize statistics for this c value
        stats_for_c = {
            'length_mismatches': 0, # Count of b values rejected due to length
            'non-monotonic': 0,     # Diff was not monotonic decreasing
            'internal_zeros': 0,    # Count of diff lists with internal zeros
            'last_exp_too_large': 0, # Count of diff lists with last exp > 1
            'b_attempts': 0,        # Count of b values that passed length check
            'a_attempts': 0,        # Number of times it trys to find a!
            'smaller_rejects': 0,   # Count of a values rejected for small exponents
            'PF-slope': 0           # Difference between exponent for 2 and 5
        }
        
        c_factors = searcher.get_factorial_factors(c_idx)
        c_length = len(c_factors) # How long is this factorization
        stats_for_c['c_list_length'] = c_length

        if not searcher.use_progress_bar and (c_idx + 1 == end_idx) and not quiet:
            print(f"Finished check for {c_idx}! with {c_length} prime factors")

        # Use the difference between the 2 and 5 exponents aa a 
        # proxy for the exponent curve slope
        if c_length > 2:
            stats_for_c['PF-slope'] = c_factors[0] - c_factors[2]

        # Calculate b range size
        first_b = c_idx - 2 # c_idx - 1 is the trivial case b_idx
        last_b = max(2, c_idx-MAX_B_LOOKBACK)
    
        
        for b_idx in range(first_b, last_b, -1):
            b_factors = searcher.get_factorial_factors(b_idx)
            stats_for_c['b_attempts'] += 1

            # If the user specified a number of lower primes to use
            # in the search, stup using haigher factors of b_factors
            not_using_k_primes = (k_primes == 0) or (len(b_factors) <= k_primes)
            # print(f"not_using_k_primes is {not_using_k_primes} and len(b_factors) is {len(b_factors)}")

            if  not_using_k_primes: # Do these tests for full search
                # Lengths must match unless we are only using the lower k primes
                if len(b_factors) != c_length: # Lengths must match
                    stats_for_c['length_mismatches'] += 1
                    break # Go get next c_idx
            
                # With matching length b_factors, we need to find the highest
                # order prime exponents that don't subtract to zero
                # and check that subtraction equals 1
                d_idx = c_length-2
                for j in range(c_length-2, -1, -1):
                    d_idx = j # Find the high prime index for diff
                    if c_factors[j] > b_factors[j]:
                        break

                # Top exp of diff must == 1
                if c_factors[d_idx] - b_factors[d_idx] > 1: 
                    stats_for_c['last_exp_too_large'] += 1
                    break # go to next c_idx 

            else: # Here if we are using k_primes part of b_factors
                d_idx = k_primes-1 # That is as far as we go into b_factors
                    
           
            diff = [] # Start building the difference in prime factor exponetns
            diff.append(c_factors[0] - b_factors[0]) # Do the 2 power first
            if diff[0] == 0: # Cannot be zero
                stats_for_c['internal_zeros'] += 1 # mark in stats
                continue # move on to next b_idx
            error_flag = False

            # if c_idx < 11:
            #    print(f"b_idx = {b_idx}, d_idx = {d_idx} and length of b_factors is {len(b_factors)}")
            #    print(f"b_factors are {b_factors} and c_factors are {c_factors}")
                

            for j in range(1,d_idx+1): # Start at power of 3  
                # print(f"j = {j}")                     
                diff.append(c_factors[j] - b_factors[j]) # Subtract power list

                if diff[j] == 0: # Cannot be zero
                    stats_for_c['internal_zeros'] += 1 # mark in stats
                    error_flag = True
                    break # monve on to the next b_idx

                if diff[j] > diff[j-1]: # Powers must be monotnic decreasing
                    stats_for_c['non-monotonic'] += 1
                    error_flag = True
                    break # move on to the next b_idx
            # if c_idx < 11:
            #    print(f"diff = {diff}")
                
            if error_flag:
                continue # Move on to next b_idx
            
 
            # At this point we have a proper string in diff
            # that could match a factorial
            # print(f"Attempting to find a! for {c_idx}! and {b_idx}! b = c - {c_idx-b_idx} and length of diff is {len(diff)}")
            # print(f"c_factors = {c_factors}")
            # print(f"b_factors = {b_factors}")
            # print(f"Diff = {diff}")
            # Search backwards from b-1 for matching a

            # Assert that window hasn't slid before searching for small a_idx value
            # print(f"Window start = {searcher.window_start}")
            assert searcher.window_start <= 2, (
                "Window has slid but attempting to search low a values. "
                "Initial window size must be at least 2200."
                )
            
            stats_for_c['a_attempts'] += 1  # Count the attempt

            # Search upward from 2 for matching a
            error_flag = False
            for a_idx in range(2, b_idx):
                a_factors = searcher.get_factorial_factors(a_idx)
    
                if len(a_factors) > len(diff) and not_using_k_primes:
                    break  # All higher a values will be too large
    
                if len(a_factors) < len(diff):
                    continue  # Try next a value

                # Check each exponent in order
                for i in range(len(diff)):
                    if a_factors[i] < diff[i]:
                        error_flag = True
                        stats_for_c['smaller_rejects'] += 1
                        break
                    if a_factors[i] > diff[i]:
                        break
                else:  # All exponents match
                    exceptions.append((a_idx, b_idx, c_idx))
                    break
                
            if error_flag:
                break
        
        # Record stats for this c value if we had any activity
        if searcher.stats and (stats_for_c['b_attempts'] > 0 or
                             stats_for_c['length_mismatches'] > 0):
            searcher.stats.record_c_stats(c_idx, stats_for_c)
    
    return exceptions

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Factorial Product Search Program - Searches for solutions to a! x b! = c!',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--limit', type=int, default=SEARCH_LIMIT,
                       help='Upper limit for factorial search')
    parser.add_argument('--window-size', type=int, default=WINDOW_SIZE,
                       help='Size of sliding window for factorial lists')
    parser.add_argument('--no-stats', action='store_true',
                       help='Disable statistics collection')
    parser.add_argument('--progress-bar', action='store_true',
                       help='Use progress bar instead of interval reporting')
    parser.add_argument('--k_primes', type=int, default=0,
                       help='Only use first k prime factors in search')
    
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_arguments()
    STATS_ENABLED = not args.no_stats
    searcher = FactorialSearch(limit=args.limit+1, use_progress_bar=args.progress_bar, k_primes=args.k_primes)
    searcher.search()

# END OF PROGRAM

