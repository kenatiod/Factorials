# Factorials
Code to explore a! x b! = c! where 1 < a < b < c-1

I became interested in the solutions to this factorial equation. It turns out that
there is a countable infinity of solutions of the form:
n! x (n!-1)! = (n!)! which comes from the definition k x (k-1)! = k! if k = n!.

However it is also known that 6! x 7! = 10! which doe not fit that
general pattern. The question is why? and are there any other
solutions that don't fit? Can we prove that there are not?

The program Factorial_Statistics.py is a command line python program that
searches the integer solution space for triplets that solve the
equation, but don't follow the usual pattern where b! = (c-1)! resulting in
a! = c. It does this by succesive tests of integers for c, and then
looping back from b = c-2 to c/2, checks to see if c!/b! equals
any a! for a < b.

Run "python Factorial_Statistics.py" to do a default search up to 2500! and 
generate a .csv statistics file on the way. The statistics are useful
in a "vanishing possibility" proof that (6, 7, 10) is the only special
solution.

Run "python Factorial_Statistics.py --limit 10000000 --no-stats" to search
for solutions up to c = 1,000,000 with no statistics file generated (because
it would be huge) or even higher if your machine will do it.

If you want to watch it make progress while it runs, use --progress-bar to
see how fast it is generating lists of primes and how fast it
is searching those lists.

The code uses the python tqdm module for the console parser. You may need to:
pip install tqdm

Using this program I was able to find a special number, where c = 2080, such
that for all c > 2080, no candidate for b < c-1 and > c/2 can generate
a c!/b! that could possibly be equal to a!

-Ken Clements

