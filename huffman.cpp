/*
 * Christopher Lee (2022_08_07)
 *
 * The file implements both the encoding and
 * decoding of the Huffman Encoding algorithm
 */

#include "bits.h"
#include "treenode.h"
#include "huffman.h"
#include "map.h"
#include "vector.h"
#include "priorityqueue.h"
#include "strlib.h"
#include "testing/SimpleTest.h"
using namespace std;

/**
 * Given a Queue<Bit> containing the compressed message bits and the encoding tree
 * used to encode those bits, decode the bits back to the original message text.
 *
 * You can assume that tree is a well-formed non-empty encoding tree and
 * messageBits queue contains a valid sequence of encoded bits.
 *
 * Your implementation may change the messageBits queue however you like. There
 * are no requirements about what it should look like after this function
 * returns. The encoding tree should be unchanged.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 *
 * @brief This function decodes a given message by checking if a given prefix cooresponds to
 * a letter in the encoding tree, which it then appends to a string
 * @return the decoded string
 *
 * @param tree the encoding tree
 * @param messageBits a sequence of Bits that represents the compressed message
 */
string decodeText(EncodingTreeNode* tree, Queue<Bit>& messageBits) {
    /* TODO: Implement this function. */
    string text = "";
    EncodingTreeNode* currentNode = tree;
    while (!messageBits.isEmpty()) {
        if (messageBits.dequeue() == 0) {
            currentNode = currentNode->zero;
        } else {
            currentNode = currentNode->one;
        }

        // If the prefix maps to a character, add the character, wipe the prefix
        if (currentNode->isLeaf()) {
            text += currentNode->getChar();
            currentNode = tree;
        }
    }
    return text;
}

/**
 * Reconstruct encoding tree from flattened form Queue<Bit> and Queue<char>.
 *
 * You can assume that the queues are well-formed and represent
 * a valid encoding tree.
 *
 * Your implementation may change the queue parameters however you like. There
 * are no requirements about what they should look like after this function
 * returns.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 *
 * @brief This function converts a sequence of Bits and Characters and converts them into
 * and encoding tree.
 *
 * @return proper Encoding Tree
 *
 * @param treeShape the Bits that represent the shape of the tree (0 - leaf, 1 - parent)
 * @param treeLeaves the characters that are included in the tree
 */
EncodingTreeNode* unflattenTree(Queue<Bit>& treeShape, Queue<char>& treeLeaves) {
    /* TODO: Implement this function. */
    if (treeShape.dequeue() == 0) {
        return new EncodingTreeNode(treeLeaves.dequeue());
    }

    return new EncodingTreeNode(unflattenTree(treeShape, treeLeaves), unflattenTree(treeShape, treeLeaves));
}

/**
 * Decompress the given EncodedData and return the original text.
 *
 * You can assume the input data is well-formed and was created by a correct
 * implementation of compress.
 *
 * Your implementation may change the data parameter however you like. There
 * are no requirements about what it should look like after this function
 * returns.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 *
 * @brief A function that given an object Encoded Data, returns the string
 * from wihtin
 *
 * @return decoded/decompressed string
 *
 * @param data and struct containing {sequence of bits defining encoding tree shape,
 *                              sequence of characters in the encoding tree,
 *                              sequence of bits represnting the compressed message}
 */
string decompress(EncodedData& data) {
    /* TODO: Implement this function. */
    EncodingTreeNode *root = unflattenTree(data.treeShape, data.treeLeaves);
    string text = decodeText(root, data.messageBits);
    deallocateTree(root);
    return text;
}

/**
 * Constructs an optimal Huffman coding tree for the given text, using
 * the algorithm described in lecture.
 *
 * Reports an error if the input text does not contain at least
 * two distinct characters.
 *
 * When assembling larger trees out of smaller ones, make sure to set the first
 * tree dequeued from the queue to be the zero subtree of the new tree and the
 * second tree as the one subtree.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 *
 * @brief This function builds a huffman tree from a given string
 *
 * @return a Huffman tree that encodes a given string "text"
 *
 * @param text the string that is to be compressed
 */
EncodingTreeNode* buildHuffmanTree(string text) {
    /* TODO: Implement this function. */
    Map<char, int> occurrences;         // Histogram
    for (char ch : text) {
        occurrences[ch]++;
    }

    if (occurrences.size() < 2) {
        return nullptr;
    }

    PriorityQueue<EncodingTreeNode*> pq;

    for (char ch : occurrences.keys()) {
        pq.enqueue(new EncodingTreeNode(ch), occurrences[ch]);
    }

    while (pq.size() > 1) {
        int priority = pq.peekPriority();
        EncodingTreeNode* left = pq.dequeue();

        priority += pq.peekPriority();
        EncodingTreeNode* right = pq.dequeue();

        pq.enqueue(new EncodingTreeNode(left, right), priority);
    }

    return pq.dequeue();    // only element in pq
}

/**
 * @brief encodingTableHelper fills out the map that relates a character to a prefix of bits
 *          by tracing all paths of the tree to find cooresponding characters
 * @param table is the map that is being filled out
 * @param tree is the root node of the Huffman tree
 * @param path is the current path the function is on
 */
void encodingTableHelper(Map<char, string>& table, EncodingTreeNode* tree, string path) {
    if (tree->isLeaf()) {
        table.put(tree->getChar(), path);
        return;
    }

    encodingTableHelper(table, tree->zero, path + "0");
    encodingTableHelper(table, tree->one, path + "1");
}
/**
 * @brief encodingTable returns a map relating a character to a prefix of bits
 * @param tree the root node of the Huffman tree
 * @return A map relating a character to a prefix of bits
 */
Map<char, string> encodingTable(EncodingTreeNode* tree) {
    Map<char, string> table;
    if (!tree) {
        return table;
    }

    encodingTableHelper(table, tree, "");
    return table;
}

/**
 * Given a string and an encoding tree, encode the text using the tree
 * and return a Queue<Bit> of the encoded bit sequence.
 *
 * You can assume tree is a valid non-empty encoding tree and contains an
 * encoding for every character in the text.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 * @brief encodeText transforms a string of text into a sequence of bits
 *          that represents the original string
 *
 * @param tree the root of the Huffman tree
 * @param text the string that is to be encoded
 *
 * @return a sequence of bits that represents the string
 */
Queue<Bit> encodeText(EncodingTreeNode* tree, string text) {
    /* TODO: Implement this function. */
    Queue<Bit> queue;
    Map<char, string> table = encodingTable(tree);
    for (char ch : text) {
        for (char bit : table[ch]) {
            if (bit == '0') {
                queue.enqueue(0);
            } else {
                queue.enqueue(1);
            }
        }
    }
    return queue;
}

/**
 * Flatten the given tree into a Queue<Bit> and Queue<char> in the manner
 * specified in the assignment writeup.
 *
 * You can assume the input queues are empty on entry to this function.
 *
 * You can assume tree is a valid well-formed encoding tree.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 *
 * @brief flattenTree converts a Huffman tree in memory to a sequence of Bit and characters
 *          that may be used to rebuild the tree (preorder traversal)
 * @param tree the root of the Huffman Tree
 * @param treeShape the sequence of Bit representing the shape of the tree (filled by the funciton)
 * @param treeLeaves the sequence of characters contains the characters in the tree (filled by the function)
 */
void flattenTree(EncodingTreeNode* tree, Queue<Bit>& treeShape, Queue<char>& treeLeaves) {
    /* TODO: Implement this function. */
    // Should only be triggered if the entire tree is empty
    if (!tree) {
        return;
    }

    if (tree->isLeaf()) {
        treeShape.enqueue(0);
        treeLeaves.enqueue(tree->getChar());
        return;
    }

    treeShape.enqueue(1);
    flattenTree(tree->zero, treeShape, treeLeaves);
    flattenTree(tree->one, treeShape, treeLeaves);
}

/**
 * Compress the input text using Huffman coding, producing as output
 * an EncodedData containing the encoded message and flattened
 * encoding tree used.
 *
 * Reports an error if the message text does not contain at least
 * two distinct characters.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 * @brief This function utilizes other functions to compress a given string into comononts of data that
 *          are meant to compress the string.
 *
 * @param messageText the string that is meant to be compressed
 *
 * @return a struct containing {sequence of bits defining encoding tree shape,
 *                              sequence of characters in the encoding tree,
 *                              sequence of bits represnting the compressed message}
 */
EncodedData compress(string messageText) {
    /* TODO: Implement this function. */
    Queue<Bit> treeShape;
    Queue<char> treeLeaves;

    EncodingTreeNode* tree = buildHuffmanTree(messageText);
    flattenTree(tree, treeShape, treeLeaves);
    EncodedData data = EncodedData{treeShape, treeLeaves, encodeText(tree, messageText)};

    deallocateTree(tree);
    return data;
}

/* * * * * * Testing Helper Functions Below This Point * * * * * */
/**
 * @brief createExampleTree creates the tree listed below
 * @return The tree listed below
 */
EncodingTreeNode* createExampleTree() {
    /* Example encoding tree used in multiple test cases:
     *                *
     *              /   \
     *             T     *
     *                  / \
     *                 *   E
     *                / \
     *               R   S
     */
    /* TODO: Implement this utility function needed for testing. */

    EncodingTreeNode* root = new EncodingTreeNode(new EncodingTreeNode('R'), new EncodingTreeNode('S'));

    root = new EncodingTreeNode(root, new EncodingTreeNode('E'));

    root = new EncodingTreeNode(new EncodingTreeNode('T'), root);

    return root;

}

/**
 * @brief deallocateTree deletes any memory that the tree uses on the heap (postorder traversal)
 * @param t the root of the encoding tree
 */
void deallocateTree(EncodingTreeNode* t) {
    /* TODO: Implement this utility function needed for testing. */
    if (!t) {
        return;
    }

    deallocateTree(t->zero);
    deallocateTree(t->one);
    delete t;
}

/**
 * @brief areEqual recurses the same path on a and b to ensure they are identical
 * @param a root of tree a
 * @param b root of tree b
 * @return whether a and b are identical to each other
 */
bool areEqual(EncodingTreeNode* a, EncodingTreeNode* b) {
    /* TODO: Implement this utility function needed for testing. */
    // Checks if they are both empty
    if (!a || !b) {
        return a == b;
    }

    // Checks leaf nodes with the same path
    if (a->isLeaf() && b->isLeaf()) {
        return a->getChar() == b->getChar();
    }

    // Trverses the Tree with the same path
    return areEqual(a->zero, b->zero) && areEqual(a->one, b->one);
}

/* * * * * * Test Cases Below This Point * * * * * */

/* TODO: Write your own student tests. */


STUDENT_TEST("Test multiples end-to-end compress -> decompress (testing if its lossless)") {
    Vector<string> inputs = {
        "HAPPY HIP HOP",
        "Nana Nana Nana Nana Nana Nana Nana Nana Batman"
        "Research is formalized curiosity. It is poking and prying with a purpose. – Zora Neale Hurston",
    };

    for (string input: inputs) {
        EncodedData one = compress(input);
        EncodedData two = compress(one.messageBits.toString());
        EncodedData three = compress(two.messageBits.toString());


        string output = decompress(three);
        EXPECT_EQUAL(two.messageBits.toString(), output);
        output = decompress(two);
        EXPECT_EQUAL(one.messageBits.toString(), output);
        output = decompress(one);
        EXPECT_EQUAL(input, output);


        // further compressions have no benefit either
        EXPECT((three.messageBits.toString() + three.treeLeaves.toString() + three.treeShape.toString()).length() >=
               (two.messageBits.toString() + two.treeLeaves.toString() + two.treeShape.toString()).length());
        EXPECT((two.messageBits.toString() + two.treeLeaves.toString() + two.treeShape.toString()).length() >=
               (one.messageBits.toString() + one.treeLeaves.toString() + one.treeShape.toString()).length());
    }
}

STUDENT_TEST("Test compression for strings with many distinct characters") {
    string input = "qwertyuiopasdfghjkzxcvbnm,[];,.4567";
    Vector<string> inputs = {
        "ABC",
        "qwertyuiopasdfghjkzxcvbnm,[];,.4567",
        "The quick brown fox jumps over the lazy dog"
    };


    for (string input: inputs) {

        EncodedData data = compress(input);
        // Every character contains 8 bits
        EXPECT((data.messageBits.toString() + data.treeLeaves.toString() + data.treeShape.toString()).length() > input.length() * 8);
    }
}




/* * * * * Provided Tests Below This Point * * * * */

PROVIDED_TEST("decodeText, small example encoding tree") {
    EncodingTreeNode* tree = createExampleTree(); // see diagram above
    EXPECT(tree != nullptr);

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(decodeText(tree, messageBits), "E");

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(decodeText(tree, messageBits), "SET");

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1}; // STREETS
    EXPECT_EQUAL(decodeText(tree, messageBits), "STREETS");

    deallocateTree(tree);
}

PROVIDED_TEST("unflattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  treeShape  = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'T', 'R', 'S', 'E' };
    EncodingTreeNode* tree = unflattenTree(treeShape, treeLeaves);

    EXPECT(areEqual(tree, reference));

    deallocateTree(tree);
    deallocateTree(reference);
}

PROVIDED_TEST("decompress, small example input") {
    EncodedData data = {
        { 1, 0, 1, 1, 0, 0, 0 }, // treeShape
        { 'T', 'R', 'S', 'E' },  // treeLeaves
        { 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 } // messageBits
    };

    EXPECT_EQUAL(decompress(data), "TRESS");
}

PROVIDED_TEST("buildHuffmanTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    EncodingTreeNode* tree = buildHuffmanTree("STREETTEST");
    EXPECT(areEqual(tree, reference));

    deallocateTree(reference);
    deallocateTree(tree);
}

PROVIDED_TEST("encodeText, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(encodeText(reference, "E"), messageBits);

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(encodeText(reference, "SET"), messageBits);

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1 }; // STREETS
    EXPECT_EQUAL(encodeText(reference, "STREETS"), messageBits);

    deallocateTree(reference);
}

PROVIDED_TEST("flattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  expectedShape  = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> expectedLeaves = { 'T', 'R', 'S', 'E' };

    Queue<Bit>  treeShape;
    Queue<char> treeLeaves;
    flattenTree(reference, treeShape, treeLeaves);

    EXPECT_EQUAL(treeShape,  expectedShape);
    EXPECT_EQUAL(treeLeaves, expectedLeaves);

    deallocateTree(reference);
}

PROVIDED_TEST("compress, small example input") {
    EncodedData data = compress("STREETTEST");
    Queue<Bit>  treeShape   = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeChars   = { 'T', 'R', 'S', 'E' };
    Queue<Bit>  messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0 };

    EXPECT_EQUAL(data.treeShape, treeShape);
    EXPECT_EQUAL(data.treeLeaves, treeChars);
    EXPECT_EQUAL(data.messageBits, messageBits);
}

PROVIDED_TEST("Test end-to-end compress -> decompress") {
    Vector<string> inputs = {
        "HAPPY HIP HOP",
        "Nana Nana Nana Nana Nana Nana Nana Nana Batman"
        "Research is formalized curiosity. It is poking and prying with a purpose. – Zora Neale Hurston",
    };

    for (string input: inputs) {
        EncodedData data = compress(input);
        string output = decompress(data);

        EXPECT_EQUAL(input, output);
    }
}
