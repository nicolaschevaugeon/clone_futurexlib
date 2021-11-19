Table of Contents
=================
* Naming convention
* Implement design
* Committing

Naming convention
=================
All function and class are supposed to use a name without '_' and starting by
a lower case letter. For every word embedded in the name the first letter is
then in upper case.
For example a function that count odd number will be named for example :
countNumber

All instance are suppose to use lower case letter with  '_' to separate
word.

For example in countNumber a local variable that is used to count information
may be call :
count_number

Type aliases (typedef or using) must end with '_t' so that we imediatelly know that it is a type alias. 
Exceptions :
* `iterator`
* `const_iterator`
* `size_type`
* all aliases already defined in the `stl` library

For example in countNumber the type of count_number may be call :
idx_t

The bunch of instruction is then
typedef long int idx_t;
idx_t countNumber(....
{
  idx_t count_number=0;

  ....

  return count_number;
}

Implement design
================

Committing
==========
Before committing run test cases (Xtest) and check that your feature does not
break anything. 

For commit always use the name of the folder in Xfiles where you have done 
modification so that any commit may immediately be identified.
For example if you change something in Octree lib you start you commit message
by "Octree:".

Avoid commit modification on many library at the same times. Split your commit 
into small commit that deal with one folder only. This offer possibility to 
revert from modifications on one library independently of the others.

Use as much as possible local branch to add your feature with as many commit 
has you want. But limit, with a git merge --squash, the number of commits
to official branch only to important step of your feature's implementation. 
Don't squeeze thing to much thought as it is always interesting to have not to 
big commit to chase buggs.