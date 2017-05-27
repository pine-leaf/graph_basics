from random import random
import os

def build_random_graph(L, p = 0.5, s = 0.0):
  degree_array = [L-1]
  edges_array = []
  weights_array = []
  for j in range(1,L):
    edges_array += [j]
    weights_array += [0.0]
  for i in range(1,L):
    degree_array += [0]
    for j in range(1,L):
      if i == j: continue
      if not random() < p: continue
      degree_array[i] += 1
      edges_array  += [j]
      weights_array += [random() - s]
  graph_file = open("Graph.in","w")
  graph_file.write("V: " + str(L) + "\n")
  graph_file.write("E: " + str(len(edges_array)) + "\n")
  next_line = "Degrees: "
  for d in degree_array:
    next_line += str(d) + " "
  graph_file.write(next_line + "\n")
  next_line = "Edges: "
  print("There are " + str(len(edges_array)) + " edges.")
  for e in edges_array:
    next_line += str(e) + " "
  graph_file.write(next_line + "\n")
  next_line = "Weights: "
  for w in weights_array:
    next_line += str(w) + " "
  graph_file.write(next_line + "\n")
  graph_file.close()

def test_params(L, p = 0.5, s = 0.0):
  build_random_graph(L, p, s)
  print("Graph built...")
  os.system('./main.o')

test_params(500, 0.1, 0.1)

def build_latin_square(D = 4):
  N = D**2
  degree_array = []
  edges_array = []
  weights_array = []
  for i in range(N):
    degree_array += [0]
    for j in range(N):
      if i == j: continue
      if (i - j) % D == 0 or int(i/D) == int(j/D):
        #print("(" + str(i) + "," + str(j) + ")")
        degree_array[i] += 1
        edges_array += [j]
        weights_array += [1.0]
  graph_file = open("Graph.in","w")
  graph_file.write("V: " + str(N) + "\n")
  graph_file.write("E: " + str(len(edges_array)) + "\n")
  next_line = "Degrees: "
  for d in degree_array:
    next_line += str(d) + " "
  graph_file.write(next_line + "\n")
  next_line = "Edges: "
  print("There are " + str(len(edges_array)) + " edges.")
  for e in edges_array:
    next_line += str(e) + " "
  graph_file.write(next_line + "\n")
  next_line = "Weights: "
  for w in weights_array:
    next_line += str(w) + " "
  graph_file.write(next_line + "\n")
  graph_file.close()


build_latin_square()

def get_latin_squares():
  inFile = open('debug.out','r')
  squares = []
  square = []
  row = []
  #current_square = []
  for L in inFile:
    L = L.rstrip()
    L = L.split()
    if not len(L): continue
    if L[0] == "ERROR:": continue
    if L[0] == "Colored:": continue
    if L[0] == "Latin":
      if square:
        squares += [square]
      square = []
      continue
    row = []
    for w in L:
      w = int(w)
      row += [w]
    square += [row]
  if square:
    squares += [square]
  return squares

S = get_latin_squares()

def count_pairs(squares):
  if not len(squares) == 2:
    return
  pairs = {}
  S1 = squares[0]
  S2 = squares[1]
  for i in range(len(S1)):
    for j in range(len(S1[0])):
      p = (S1[i][j],S2[i][j])
      if not p in pairs:
        pairs[p] = 1
      else:
        pairs[p] += 1
        print(str(i) + ", " + str(j) + ": " + str(p) + " -> " + str(pairs[p]))
  #for p in pairs.keys():
  #  print(str(p) + ": " + str(pairs[p]))
  return pairs

count_pairs(S)

def xor_latin_square(N = 4):
  file = open('latin_square.in','w')
  file.write("D: " + str(2**N) + "\n")
  for i in range(2**N):
    line = ""
    for j in range(2**N):
      s = (i^j) + 1
      line += str(s) + " "
    file.write(line+"\n")
  file.close()

xor_latin_square()

def sum_latin_square(N = 2):
  file = open('latin_square.in','w')
  file.write("D: " + str(2*N+1) + "\n")
  for i in range(2*N+1):
    line = ""
    for j in range(2*N+1):
      s = (i+j)%(2*N+1) + 1
      line += str(s) + " "
    file.write(line + "\n")
  file.close()

sum_latin_square()


