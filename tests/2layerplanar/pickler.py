import pickle

with open('test.pickle', 'rb') as f:
  dp = pickle.load(f)

print(dp)
