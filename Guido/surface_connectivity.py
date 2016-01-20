#!/usr/bin/env python
__author__ = 'Guido Falk von Rudorff'
import sys

if len(sys.argv) != 2:
	raise ValueError('Usage: script surface_patterns_file')

pattern_file = sys.argv[1]
transitions = dict()
transitions['start'] = dict()
transitions['stop'] = dict()
starting = ['start',]
stopping = ['stop',]
with open(pattern_file) as fh:
	lastkey = ''
	lastpattern = ''
	for line in fh:
		traj, surface, frame, pattern = line.strip().split()
		thiskey = traj + surface
		if lastpattern != '' and lastpattern not in transitions:
			transitions[lastpattern] = dict()
		if thiskey != lastkey:
			lastkey = thiskey
			if lastpattern != '':
				transitions[lastpattern]['stop'] = 1
				stopping.append(lastpattern)
			lastpattern = pattern
			starting.append(pattern)
			transitions['start'][pattern] = 1
			continue

		if pattern not in transitions[lastpattern]:
			transitions[lastpattern][pattern] = 0
		transitions[lastpattern][pattern] += 1
		lastpattern = pattern

extras = dict()
extras['start'] = 'shape=box fillcolor=green style=filled'
extras['stop'] = 'shape=box fillcolor=red style=filled'
for origin in transitions:
	if origin not in extras:
		extras[origin] = ''
	if origin in transitions[origin] and transitions[origin][origin] > 1000:
		extras[origin] += ' fillcolor=orange style=filled'

print 'digraph G {\n\trankdir=TB;\n\tsubgraph{'
for origin in transitions:
	if origin in starting or origin in stopping or origin.startswith('s'):
		continue
	print '\t\t%s [label="%s" %s]' % (origin, origin, extras[origin])
print '\t\t {\n\t\t\trank=source;',
for origin in starting:
	print '%s [%s];' % (origin, extras[origin]),
print '\t\t}'
print '\t\t{rank=sink;',
for origin in stopping:
	print '\t\t\t%s [%s];' % (origin, extras[origin]),
print '\t\t}'
#print '\t\t{rank=source;',
#print '%s [%s]' % ('start', extras['start']),
#print '}'
#print '\t\t{rank=sink;',
#print '%s [%s]' % ('stop', extras['stop']),
#print '}'

for origin in transitions:
	if origin == 'start':
		continue
	for target in transitions[origin]:
		if target == 'stop':
			continue
		extra = ''
		print '\t\t%s -> %s [label="%d" %s]' % (origin, target, transitions[origin][target], extra)
print 'start -> stop'
print '\t}\n}'