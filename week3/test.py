import numpy as np
import time

# Conditional imports based on environment
__codon__ = False  # Toggle this for Codon testing

if __codon__:
    from biotite_codon import upgma, neighbor_joining
else:
    import biotite.sequence.phylo as phylo
    upgma = phylo.upgma
    neighbor_joining = phylo.neighbor_joining

# Define @test decorator for compatibility
def test(func):
    """Test decorator for compatibility with both Python and Codon"""
    return func

@test
def test_distances():
    """Test phylogenetic distance calculations"""
    start_time = time.time()
    
    # Create a sample distance matrix for UPGMA
    dist_matrix = np.array([
        [0.0, 0.1, 0.2, 0.3],
        [0.1, 0.0, 0.15, 0.25],
        [0.2, 0.15, 0.0, 0.35],
        [0.3, 0.25, 0.35, 0.0]
    ])
    
    # Create tree via UPGMA
    tree = upgma(dist_matrix)
    
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    if len(tree.leaves) > 1:
        dist = tree.root.distance_to(tree.leaves[0])
        for leaf in tree.leaves:
            assert leaf.distance_to(tree.root) == dist
    
    end_time = time.time()
    runtime_ms = (end_time - start_time) * 1000
    print(f"test_distances: {runtime_ms:.2f} ms")
    print("  ✓ distance calculations")

@test
def test_upgma():
    """Test UPGMA clustering algorithm"""
    start_time = time.time()
    
    # Sample distance matrix
    dist_matrix = np.array([
        [0.0, 0.1, 0.2, 0.3],
        [0.1, 0.0, 0.15, 0.25],
        [0.2, 0.15, 0.0, 0.35],
        [0.3, 0.25, 0.35, 0.0]
    ])
    
    # Create tree using UPGMA
    tree = upgma(dist_matrix)
    
    # Basic validation - tree should have correct number of leaves
    assert len(tree.leaves) == len(dist_matrix)
    
    # Validate tree structure
    assert tree.root is not None
    
    end_time = time.time()
    runtime_ms = (end_time - start_time) * 1000
    print(f"test_upgma: {runtime_ms:.2f} ms")
    print("  ✓ UPGMA algorithm")

@test
def test_neighbor_joining():
    """Test neighbor-joining tree construction"""
    start_time = time.time()
    
    # Distance matrix from the original test
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ])
    
    # Create tree using neighbor-joining
    # Execute all tests
test_distances()
test_upgma()
test_neighbor_joining()

if __name__ == "__main__":
    print("=== Phylogenetic Algorithm Tests ===")
    env = "Codon" if __codon__ else "Python"
    print(f"Testing {env} implementation")
    print()
    
    overall_start = time.time()
    
    tests = [test_distances, test_upgma, test_neighbor_joining]
    test_results = []
    
    for test_func in tests:
        try:
            test_func()
            test_results.append((test_func.__name__, "PASSED", None))
        except Exception as e:
            test_results.append((test_func.__name__, "FAILED", str(e)))
            print(f"❌ {test_func.__name__} FAILED: {e}")
    
    overall_end = time.time()
    total_runtime_ms = (overall_end - overall_start) * 1000
    
    print()
    print("=== Test Results Summary ===")
    
    passed_count = 0
    for test_name, status, error in test_results:
        if status == "PASSED":
            print(f"✅ {test_name}: {status}")
            passed_count += 1
        else:
            print(f"❌ {test_name}: {status} - {error}")
    
    print(f"\nTests passed: {passed_count}/{len(tests)}")
    print(f"Total runtime: {total_runtime_ms:.2f} ms")
    
    # Exit with error code if any tests failed
    if passed_count < len(tests):
        exit(1)
    else:
        print("🎉 All tests passed!")
        exit(0)
    test_tree = phylo.neighbor_joining(dist)
    
    # Basic validation - tree should have correct number of leaves
    assert len(test_tree.leaves) == len(dist)
    
    # Validate tree structure
    assert test_tree.root is not None
    
    end_time = time.time()
    runtime_ms = (end_time - start_time) * 1000
    print(f"test_neighbor_joining: {runtime_ms:.2f} ms")
    print("  ✓ Neighbor Joining algorithm")

if __name__ == "__main__":
    print("=== Phylogenetic Algorithm Tests ===")
    print("Testing Biotite phylo package functionality")
    print()
    
    overall_start = time.time()
    
    tests = [test_distances, test_upgma, test_neighbor_joining]
    test_results = []
    
    for test_func in tests:
        try:
            test_func()
            test_results.append((test_func.__name__, "PASSED", None))
        except Exception as e:
            test_results.append((test_func.__name__, "FAILED", str(e)))
            print(f"❌ {test_func.__name__} FAILED: {e}")
    
    overall_end = time.time()
    total_runtime_ms = (overall_end - overall_start) * 1000
    
    print()
    print("=== Test Results Summary ===")
    
    passed_count = 0
    for test_name, status, error in test_results:
        if status == "PASSED":
            print(f"✅ {test_name}: {status}")
            passed_count += 1
        else:
            print(f"❌ {test_name}: {status} - {error}")
    
    print(f"\nTests passed: {passed_count}/{len(tests)}")
    print(f"Total runtime: {total_runtime_ms:.2f} ms")
    
    # Exit with error code if any tests failed
    if passed_count < len(tests):
        exit(1)
    else:
        print("🎉 All tests passed!")
        exit(0)
