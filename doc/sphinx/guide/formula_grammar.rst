.. _formula:

***************
Formula Grammar
***************

Chemical formula strings can be quite rich.  While they do not provide
any information about molecular structure, they do all complete control
over chemical composition. 

* Formula strings consist of counts and atoms such as "CaCO3+6H2O".          

* Groups can be separated by '+' or space, so "CaCO3 + 6H2O" works as well. 

* Groups can be defined using parentheses, such as "CaCO3 (H2O)6". 

* Parentheses can nest: "(CaCO3(H2O)6)1". 

* Isotopes are represented by index, e.g., "CaCO[18]3+6H2O". 

* Counts can be integer or decimal, e.g. "CaCO3+(3HO1.5)2".


Formulas are parsed from strings using the following grammar::

        number    :: [1-9][0-9]*
        fraction  :: ( | '0' | number) '.' [0-9]*
        count     :: number | fraction | ''
        symbol    :: [A-Z][a-z]*
        isotope   :: '[' number ']' | ''
        element   :: symbol isotope count
        separator :: '+' | ' '
        group     :: count element+ | '(' formula ')' count
        grammar   :: group separator formula | group | ''

