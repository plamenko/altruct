#pragma once

#include "algorithm/collections/collections.h"

#include <queue>
#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * A trie structure post-processed with the Aho–Corasick algorithm.
 *
 * Note, if the ALPHABET_SIZE is big (e.g. 256 letters, 8 bit),
 * one can split each letter into two or more symbols (e.g. 4 bits each).
 *
 * Space complexity: `O(ALPHABET_SIZE * DICT_SIZE)`.
 * Time complexities:
 *   insert:    `O(ALPHABET_SIZE * WORD_SIZE)`
 *   build:     `O(ALPHABET_SIZE * DICT_SIZE)`
 *   move_next: `O(1)`
 */
template<int ALPHABET_SIZE = 26, typename INDEX_T = int, typename ALPHA_T = char>
class aho_corasick_trie {
public:
	typedef INDEX_T index_t;
	typedef ALPHA_T alpha_t;

	static const index_t NIL = 0;
	static const index_t ROOT = 1;
	static const index_t RESERVED = 2;

	struct node_t {
		alpha_t let;
		index_t parent;
		index_t suff_link;
		index_t word_cnt;
		index_t next[ALPHABET_SIZE];
	};

	std::vector<node_t> trie;

	aho_corasick_trie() {
		trie.resize(RESERVED);
		trie[ROOT].suff_link = ROOT;
	}

	// Insert pattern word to the trie dictionary.
	template<typename It>
	index_t insert(It begin, It end) {
		typedef typename std::iterator_traits<It>::value_type T;
		return insert(begin, end, [](T let){ return (alpha_t)let; });
	}

	// Insert pattern word to the trie dictionary.
	// @param ordinal - functor that maps element to an alphabet ordinal
	//                  E.g. [](char c){ return c - 'a'; }
	template<typename It, typename F>
	index_t insert(It begin, It end, F ordinal) {
		index_t len = (index_t)std::distance(begin, end);
		collections::reserve_more(trie, len);
		index_t cur_node = ROOT;
		for (auto it = begin; it != end; ++it) {
			cur_node = _insert(cur_node, ordinal(*it));
		}
		trie[cur_node].word_cnt++;
		return cur_node;
	}

	index_t _insert(index_t cur_node, alpha_t let) {
		index_t next = trie[cur_node].next[let];
		if (next != NIL) return next;
		next = (index_t)trie.size();
		trie.push_back({ let, cur_node, NIL, 0 });
		return trie[cur_node].next[let] = next;
	}

	// Postprocesses the trie dictionary with the Aho-Corasick algorithm.
	void build() {
		std::queue<index_t> q;
		q.push((index_t)ROOT);
		while (!q.empty()) {
			index_t cur_node = q.front(); q.pop();
			index_t* next = trie[cur_node].next;
			for (alpha_t let = 0; let < ALPHABET_SIZE; let++) {
				if (next[let]) q.push(next[let]);
			}
			index_t suff_link = _get_suff_link(cur_node);
			trie[cur_node].suff_link = suff_link;
			trie[cur_node].word_cnt += trie[suff_link].word_cnt;
		}
	}

	// Counts all occurences of all dictionary words within the given string.
	template<typename It>
	int64_t count_matches(It begin, It end) {
		typedef typename std::iterator_traits<It>::value_type T;
		return count_matches(begin, end, [](T let){ return (alpha_t)let; });
	}
	
	// Counts all occurences of all dictionary words within the given string.
	// @param ordinal - functor that maps element to an alphabet ordinal
	//                  E.g. [](char c){ return c - 'a'; }
	template<typename It, typename F>
	int64_t count_matches(It begin, It end, F ordinal) {
		int64_t r = 0;
		int cur_node = ROOT;
		for (auto it = begin; it != end; ++it) {
			cur_node = move_next(cur_node, ordinal(*it));
			r += trie[cur_node].word_cnt;
		}
		return r;
	}

	// Moves to the next state given the current state and the transition letter.
	index_t move_next(index_t cur_node, alpha_t let) {
		node_t& t = trie[cur_node];
		if (t.next[let] != NIL) return t.next[let];
		if (cur_node == ROOT) return t.next[let] = ROOT;
		index_t suff_link = _get_suff_link(cur_node);
		return t.next[let] = move_next(suff_link, let);
	}

	index_t _get_suff_link(index_t cur_node) {
		node_t& t = trie[cur_node];
		if (t.suff_link != NIL) return t.suff_link;
		if (t.parent == ROOT) return t.suff_link = ROOT;
		index_t suff_link = trie[t.parent].suff_link;
		return t.suff_link = move_next(suff_link, t.let);
	}

	// Returns the index of the root node.
	index_t root() const {
		return ROOT;
	}

	// Returns the number of nodes.
	index_t size() const {
		return (index_t)trie.size();
	}

	// Accesses node by its index.
	node_t& operator[] (index_t index) {
		return trie[index];
	}
};

} // container
} // altruct
