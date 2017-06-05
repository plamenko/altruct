#pragma once

#include "altruct/algorithm/collections/collections.h"

#include <stdint.h>
#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * A tree that enumerates palindrome substrings.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build: `O(n)`
 */
template<int ALPHABET_SIZE = 26, typename INDEX_T = int, typename ALPHA_T = uint8_t>
class palindrome_tree {
	typedef INDEX_T index_t;
	typedef ALPHA_T alpha_t;
public:
	// node that represents a palindromic substring
	struct node_t {
		index_t len;      // length of this palindromic substring
		index_t pos;      // position of the first occurence of this palindromic substring within the string
		index_t cnt;      // multiplicity of this palindromic substring
		index_t depth;    // depth in the suffix chain of this node
		index_t suff;     // node-index of the largest palindromic suffix of this node
		index_t next[ALPHABET_SIZE]; // "A".next['x'] --> "xAx"
	};

	static const index_t NIL = 0;
	static const index_t NEGAT = 1;
	static const index_t EMPTY = 2;
	static const index_t RESERVED = 3;

	std::vector<alpha_t> _string;   // string of letter ordinals (e.g. 'a' is 0)
	std::vector<node_t> _nodes;     // node container
	index_t _suff;                  // node-index of the current longest palindromic suffix
	int64_t _total;                 // total number of palindromic substrings (counting multiplicities), can be quadratic in string length

	palindrome_tree() {
		_init();
	}

	template<typename It, typename F>
	palindrome_tree(It begin, It end, F ordinal) {
		_init();
		add_all(begin, end, ordinal);
	}

	void _init() {
		_string.clear();
		_nodes.resize(RESERVED);
		_nodes[NEGAT].len = -1; _nodes[NEGAT].suff = NEGAT;
		_nodes[EMPTY].len = 0; _nodes[EMPTY].suff = NEGAT;
		_suff = EMPTY;
		_total = 0;
	}

	template<typename It, typename F>
	index_t add_all(It begin, It end, F ordinal) {
		auto len = distance(begin, end);
		collections::reserve_more(_string, len);
		collections::reserve_more(_nodes, len);
		index_t c = 0;
		for (auto it = begin; it != end; ++it) {
			c += add(*it, ordinal);
		}
		return c;
	}

	template<typename T, typename F>
	index_t add(const T& t, F ordinal) {
		return add(ordinal(t));
	}
	
	index_t add(alpha_t let) {
		_string.push_back(let);
		index_t i = _find_suffix(_suff, let);
		_suff = _nodes[i].next[let];
		if (_suff != NIL) {
			_nodes[_suff].cnt++;
			_total += _nodes[_suff].depth;
			return 0;
		}
		index_t suff2 = _find_suffix2(i, let);
		_suff = (index_t)_nodes.size();
		_nodes.push_back({});
		_nodes[_suff].len = _nodes[i].len + 2;
		_nodes[_suff].pos = (index_t)_string.size() - _nodes[_suff].len;
		_nodes[_suff].cnt = 1;
		_nodes[_suff].suff = suff2;
		_nodes[_suff].depth = _nodes[suff2].depth + 1;
		_nodes[i].next[let] = _suff;
		_total += _nodes[_suff].depth;
		return 1;
	}

	index_t _find_suffix2(index_t i, alpha_t let) {
		if (i == NEGAT) return EMPTY;
		i = _find_suffix(_nodes[i].suff, let);
		return _nodes[i].next[let];
	}

	index_t _find_suffix(index_t i, alpha_t let) {
		index_t sz = (index_t)_string.size();
		while (sz < _nodes[i].len + 2 || _string[sz - _nodes[i].len - 2] != let) {
			i = _nodes[i].suff;
		}
		return i;
	}

	// This should be called only once after all elements are added!
	void propagate() {
		for (index_t i = (index_t)_nodes.size() - 1; i >= RESERVED; i--) {
			index_t suff = _nodes[i].suff;
			_nodes[suff].cnt += _nodes[i].cnt;
		}
	}

	// Returns the number of total palindromic substrings, counting their multiplicities.
	int64_t total() const {
		return _total;
	}

	// Returns the number of distinct palindromic substrings, each counted only once.
	index_t distinct() const {
		return (index_t)_nodes.size() - RESERVED;
	}

	// Returns the index of the node representing the longest palindromic suffix.
	index_t longest_suffix() const {
		return _suff;
	}

	// Returns the index of the first node.
	index_t first() const {
		return RESERVED;
	}

	// Returns the number of nodes.
	index_t size() const {
		return (index_t)_nodes.size();
	}

	// Accesses node by its index.
	node_t& operator[] (index_t index) {
		return _nodes[index];
	}

	static alpha_t ordinal_digit(char c) { return alpha_t(c - '0'); }
	static alpha_t ordinal_lower_alpha(char c) { return alpha_t(c - 'a'); }
	static alpha_t ordinal_upper_alpha(char c) { return alpha_t(c - 'A'); }
};

} // container
} // altruct
