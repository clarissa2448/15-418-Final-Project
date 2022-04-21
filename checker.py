#!/usr/bin/env python
import logging
import sys

FORMAT = '%(asctime)-15s - %(levelname)s - %(module)10s:%(lineno)-5d - %(message)s'
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=FORMAT)
LOG = logging.getLogger(__name__)

help_message = '''
usage: checker.py [-h] [-i INPUT] [-n n] [-o OUTPUT] [-e e]

Tests the correctness of the wire routing output and cost matrix

Arguments:
  -h, --help            Show this help message and exit
  -i INPUT              Input adjacency List
  -o OUTPUT             Outputted independent set
  -n n                  2^n nodes
  -E e                  e edges
'''

def parse_args():
    args = sys.argv
    if '-h' in args or '--help' in args:
        print(help_message)
        sys.exit(1)
    if '-o' not in args or '-n' not in args or '-E' not in args or '-i' not in args:
        print(help_message)
        sys.exit(1)
    parsed = {}
    parsed['output'] = args[args.index('-o') + 1]
    parsed['n'] = args[args.index('-n') + 1]
    parsed['E'] = args[args.index('-E') + 1]
    parsed['input'] = args[args.index('-i') + 1]
    return parsed


def main(args):
    val = validate(args)
    print ("Correctness: " + str(val))

def validate(args):
    # Parse Files to Construct graph and independent set
    n = int(args['n'])
    N = 2 ** n
    
    E = args['E']
    output = open(args['output'], 'r')
    lines = output.readlines()
    if len(lines) < 1:
        LOG.error('''Input file contains has less than 2 lines,
        please check for the input format''')
        return False
    independent_set = set()
    for line in lines:
        independent_set.add(int(line))
    
    input = open(args['input'], 'r')
    adj_list = []
    lines = input.readlines()
    for line in lines:
        arr = line.split()
        arr_list = set()
        for a in arr:
            arr_list.add(int(a))
        adj_list.append(arr_list)
    # Check is independent set: no two nodes are neighbors
    for u in independent_set:
        for v in independent_set:
            if u != v and v in adj_list[u]:
                LOG.error('''Not an Independent Set''')
                return False
    
    # Check that independent set is maximal
    for u in range(N):
        if u not in independent_set:
            can_add = True
            for v in adj_list[u]:
                if v in independent_set:
                    can_add = False
            if can_add:
                LOG.error('''Independent Set not Maximal''')
                return False
    return True
    

    



if __name__ == '__main__':
    main(parse_args())
