#!/usr/bin/python3
# ==============================================================================
# OpenACC Variable Declaration Consistency Checker
# ==============================================================================
# Purpose:
#   Verifies that global variable declarations in 'params.f90' match exactly
#   with OpenACC directives:
#     1. '!$ACC declare create' in 'params.f90'
#     2. '!$ACC update device' in 'setflac.f90'
#
# Usage:
#   Run this script from the 'src/' directory:
#     python3 ../util/check-acc-decl.py
#
# Parsing Requirements & Constraints:
#   - Identical Line-by-Line Layout: The script parses lines verbatim. All three
#     locations must group variables on the same lines in the same order.
#   - No Inline Comments: Comments (lines starting with '!') must not be inserted
#     in the middle of multi-line declarations/directives as it halts parsing loops.
#   - Line Continuation: All continuation lines in directives must end with '&'.
#   - Final Parenthesis Constraint: Array variables (e.g., 'dumC(4)') should not
#     be placed at the very end of a block/clause because the parser strips
#     closing parentheses differently for F90 vs. ACC directives. Always place
#     a scalar variable (e.g., 'xmodalcpx') at the end of the lists.
#
# Output:
#   - '# of lines: A B C': Counts of parsed declaration, declare, and update lines.
#   - 'Line X mismatch in acc create/update': Identifies mismatching lines.
#
# How to Resolve Mismatches:
#   1. Align variables in 'params.f90' declarations, 'params.f90' ACC declare,
#      and 'setflac.f90' ACC update device blocks.
#   2. Match line breaks and variable ordering exactly.
#   3. Ensure continuation '&' is present at the end of every non-final line.
# ==============================================================================

decli, decld, accdecl, accupd = [], [], [], []
acc_hdr = '!$ACC '

s = open('params.f90').readlines()
nline = len(s)
n = 0
while n < nline:
    line = s[n].strip()

    # declaration for integer
    head = 'integer :: '
    if line.startswith(head):
        line = line[len(head):].strip()
        while line.endswith('&'):
            decli.append(line)
            n += 1
            line = s[n].strip()
        decli.append(line) # last line of integer declaration
        #print(line, decli)

    # declaration for double
    head = 'real*8 :: '
    if line.startswith(head):
        line = line[len(head):].strip()
        while line.endswith('&'):
            decld.append(line)
            n += 1
            line = s[n].strip()
        decld.append(line) # last line of double declaration
        #print(line, decld)


    # acc declaration for integer and double
    head = acc_hdr +'declare create('
    if line.startswith(head):
        line = line[len(head):].strip()
        while line.endswith('&'): # multi-lines of acc declare
            accdecl.append(line)
            n += 1
            line = s[n][len(acc_hdr):].strip()
        accdecl.append(line.strip(')')) # last line of acc declare
        #print(line, accdecl)

    n += 1


s = open('setflac.f90').readlines()
nline = len(s)
n = 0
while n < nline:
    line = s[n].strip()

    # acc update for integer and double
    head = acc_hdr +'update device('
    if line.startswith(head):
        line = line[len(head):].strip()
        if not line.endswith('&'): # single line of acc, skip it
            n += 1
            line = s[n].strip()
            continue
        while line.endswith('&'): # multi-lines of acc update
            accupd.append(line)
            n += 1
            line = s[n][len(acc_hdr):].strip()
        accupd.append(line[:-len(') async(1)')])  # last line of acc update
    n += 1
#print(accupd)


# diff'ing
decl = decli + decld
print('Comparing the declaration of F90, ACC-declare, and ACC-update:')
print('# of lines:', len(decl), len(accdecl), len(accupd))

ndiff = 0
for n, line in enumerate(decl):
    if line != accdecl[n]:
        ndiff += 1
        print('Line', n, 'mismatch in acc create:', line, accdecl[n])

    if line != accupd[n]:
        ndiff += 1
        print('Line', n, 'mismatch in acc update:', line, accupd[n])

print('Comparison finished,', ndiff, 'mismatch(es) found.')