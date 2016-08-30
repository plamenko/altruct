#pragma once

#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * Range minimum query structure.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build: `O(n)`
 *   query: `O(1)`
 */
template<typename VALUE_T, typename COMPARE_T = std::less<VALUE_T>, typename INDEX_T = int, int BLOCK_SIZE = 8>
class direct_rmq {
public:
	typedef VALUE_T value_t;
	typedef COMPARE_T compare_t;
	typedef INDEX_T index_t;
	typedef char small_index_t;
	
	struct tree_t {
		small_index_t idx[BLOCK_SIZE][BLOCK_SIZE];
	};
	
	direct_rmq(compare_t comp = compare_t()) : comp(comp) {
		calc_ballot_numbers();
	}

	template<typename It>
	void build(It begin, It end) {
		index_t len = (index_t)std::distance(begin, end);
		blocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
		_array.assign(begin, end);
		_array.resize(blocks * BLOCK_SIZE);
		build_block_trees();
		build_sparse_table();
	}

	value_t get_value(index_t begin, index_t end) const {
		return _array.at(get_index(begin, end));
	}

	index_t get_index(index_t begin, index_t end) const {
		index_t l = begin, r = end - 1;
		index_t x = l / BLOCK_SIZE, y = r / BLOCK_SIZE, z = y - x;
		if (z == 0) return x * BLOCK_SIZE + _block_tree[x]->idx[l % BLOCK_SIZE][r % BLOCK_SIZE];
		if (z == 1) return min_index(
			x * BLOCK_SIZE + _block_tree[x]->idx[l % BLOCK_SIZE][BLOCK_SIZE - 1],
			y * BLOCK_SIZE + _block_tree[y]->idx[0][r % BLOCK_SIZE]);
		int k = log2_32(z - 2);
		return min_index(
			min_index(
				x * BLOCK_SIZE + _block_tree[x]->idx[l % BLOCK_SIZE][BLOCK_SIZE - 1],
				_sparse_table[x + 1 + blocks * k]
			),
			min_index(
				_sparse_table[y + blocks * k - (1 << k)],
				y * BLOCK_SIZE + _block_tree[y]->idx[0][r % BLOCK_SIZE]
			));
	}

	index_t size() {
		return (index_t)_array.size();
	}

private:
	compare_t comp;
	index_t blocks;
	std::vector<value_t> _array;
	std::vector<tree_t> _trees;
	std::vector<tree_t*> _block_tree;
	std::vector<index_t> _sparse_table;
	int _ballot[BLOCK_SIZE + 1][BLOCK_SIZE + 1];

	void build_block_trees() {
		value_t rp[BLOCK_SIZE + 1];
		_trees.resize(_ballot[BLOCK_SIZE][BLOCK_SIZE]);
		for (auto& t : _trees) {
			t.idx[0][0] = -1;
		}
		_block_tree.resize(blocks);
		for (index_t i = 0; i < blocks; i++) {
			// determine tree type
			const value_t *p = _array.data() + i * BLOCK_SIZE;
			int q = BLOCK_SIZE, N = 0;
			for (int j = 0; j < BLOCK_SIZE; j++) {
				while (q + j - BLOCK_SIZE > 0 && comp(p[j], rp[q + j - BLOCK_SIZE])) {
					N += _ballot[BLOCK_SIZE - j - 1][q];
					q--;
				}
				rp[q + j - BLOCK_SIZE + 1] = p[j];
			}
			_block_tree[i] = &_trees[N];
			// populate tree
			tree_t* table = _block_tree[i];
			if (table->idx[0][0] != -1) continue;
			for (int left = 0; left < BLOCK_SIZE; left++) {
				int min_idx = left;
				for (int right = left; right < BLOCK_SIZE; right++) {
					if (comp(p[right], p[min_idx])) {
						min_idx = right;
					}
					table->idx[left][right] = (small_index_t)min_idx;
				}
			}
		}
	}

	void build_sparse_table() {
		int height = 0; while (1 << height < blocks) ++height;
		_sparse_table.resize(blocks * height);
		if (height == 0) return;
		index_t *b = _sparse_table.data();
		for (index_t i = 0; i < blocks; i++) {
			b[i] = i * BLOCK_SIZE + _block_tree[i]->idx[0][BLOCK_SIZE - 1];
		}
		for (index_t t = 1; t * 2 < blocks; t *= 2) {
			std::copy(b, b + blocks, b + blocks);
			b += blocks;
			for (index_t i = 0; i < blocks - t; ++i) {
				b[i] = min_index(b[i], b[i + t]);
			}
		}
	}

	index_t min_index(index_t x, index_t y) const {
		return comp(_array[y], _array[x]) ? y : x;
	}

	void calc_ballot_numbers() {
		for (int q = 0; q <= BLOCK_SIZE; q++) {
			_ballot[0][q] = 1;
		}
		for (int p = 1; p <= BLOCK_SIZE; p++) {
			for (int q = 0; q < p; q++) {
				_ballot[p][q] = 0;
			}
			for (int q = p; q <= BLOCK_SIZE; q++) {
				_ballot[p][q] = _ballot[p][q - 1] + _ballot[p - 1][q];
			}
		}
	}

	static int log2_32(uint32_t val) {
		static const int8_t tbl[32] = {
			 0,  9,  1, 10, 13, 21,  2, 29,
			11, 14, 16, 18, 22, 25,  3, 30,
			 8, 12, 20, 28, 15, 17, 24,  7,
			19, 27, 23,  6, 26,  5,  4, 31 };
		val |= val >> 1;
		val |= val >> 2;
		val |= val >> 4;
		val |= val >> 8;
		val |= val >> 16;
		return tbl[(val * (uint32_t)0x07C4ACDD) >> 27];
	}
};

} // container
} // altruct
