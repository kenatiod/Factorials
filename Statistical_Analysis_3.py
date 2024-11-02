"""






"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def analyze_factorial_stats(stats_file):
    """Analyze factorial search statistics"""
    print("\nFactorial Search Statistics Analysis")
    print("===================================")
    
    # Read the CSV file
    df = pd.read_csv(stats_file)
    
    # Create figure for multiple plots
    plt.figure(figsize=(20, 15))
    
    # 1. Prime Factor Slope Analysis
    plt.subplot(2, 2, 1)
    plt.plot(df['c'], df['PF-slope'], 'b-', label='PF-slope')
    plt.title('Prime Factor Slope vs c')
    plt.xlabel('c')
    plt.ylabel('Slope (difference between exp of 2 and 5)')
    plt.grid(True)
    
    # 2. Attempts Analysis
    plt.subplot(2, 2, 2)
    plt.plot(df['c'], df['b_attempts'], 'r-', label='b attempts')
    plt.plot(df['c'], df['a_attempts'], 'g-', label='a attempts')
    plt.title('Search Attempts vs c')
    plt.xlabel('c')
    plt.ylabel('Number of attempts')
    plt.legend()
    plt.grid(True)
    
    # 3. Rejection Reasons
    plt.subplot(2, 2, 3)
    rejection_columns = ['length_mismatches', 'non-monotonic', 
                        'internal_zeros', 'last_exp_too_large']
    df[rejection_columns].plot(kind='area', stacked=True)
    plt.title('Rejection Reasons Over c')
    plt.xlabel('c')
    plt.ylabel('Count')
    plt.grid(True)
    
    # 4. List Length and PF-slope
    plt.subplot(2, 2, 4)
    plt.plot(df['c'], df['c_list_length'], 'r-', label='List Length')
    plt.plot(df['c'], df['PF-slope'], 'b-', label='PF-slope')
    plt.title('List Length and PF-slope vs c')
    plt.xlabel('c')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('factorial_analysis.png')
    
    # Print key statistics
    print("\nKey Statistics:")
    print("---------------")
    
    # Find last occurrences of attempts
    last_attempts = {
        'b_attempts': df[df['b_attempts'] > 0]['c'].max(),
        'a_attempts': df[df['a_attempts'] > 0]['c'].max()
    }
    print("\nLast Occurrences:")
    for attempt_type, value in last_attempts.items():
        print(f"{attempt_type}: {value}!")
    
    # Analyze PF-slope progression
    print("\nPF-slope Analysis:")
    print(f"Initial slope: {df['PF-slope'].iloc[0]}")
    print(f"Final slope: {df['PF-slope'].iloc[-1]}")
    print(f"Maximum slope: {df['PF-slope'].max()}")
    
    # Analyze list length progression
    print("\nList Length Analysis:")
    print(f"Initial length: {df['c_list_length'].iloc[0]}")
    print(f"Final length: {df['c_list_length'].iloc[-1]}")
    print(f"Maximum length: {df['c_list_length'].max()}")
    
    # Analyze rejection patterns
    print("\nRejection Pattern Analysis:")
    windows = df.groupby(df['c'] // 100).agg({
        'length_mismatches': 'sum',
        'non-monotonic': 'sum',
        'internal_zeros': 'sum',
        'last_exp_too_large': 'sum',
        'b_attempts': 'sum',
        'a_attempts': 'sum',
        'PF-slope': 'mean',
        'c_list_length': 'mean'
    })
    print("\nRejections by hundreds:")
    print(windows)

if __name__ == "__main__":
    stats_file = "factorial_search_stats_20241102_075744.csv"
    analyze_factorial_stats(stats_file)

# END OF PROGRAM
