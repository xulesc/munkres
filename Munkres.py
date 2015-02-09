import munkres
import numpy as np

class Munkres:
  def run_munkres(self, cost_matrix, do_max = 0):
    orig_shape = cost_matrix.shape
    if orig_shape[0] != orig_shape[1]:
      self.cost_matrix = np.copy(cost_matrix)
      max_dim = max(orig_shape)
      self.cost_matrix.resize(max_dim, max_dim)
    else:
      self.cost_matrix = cost_matrix
    
    c = np.array(munkres.run_munkres(self.cost_matrix, do_max)).reshape(self.cost_matrix.shape)
    return c[0:orig_shape[0], 0:orig_shape[1]]
    
if __name__ == '__main__':
  m = Munkres()
  print 'Test 1 - forces matrix copy'
  a = np.array([[1, 2, 3, 6], [2, 4, 6, 5], [3, 6, 9, 1]]).astype(np.int32)  
  print 'a = %s' %a
  print 'c = %s' %(m.run_munkres(a))
  print 'Test 2 - no matrix copy'
  a = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]]).astype(np.int32)  
  print 'a = %s' %a
  print 'c = %s' %(m.run_munkres(a))
