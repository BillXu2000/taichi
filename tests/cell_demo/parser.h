#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <stack>
#include <set>
#include <cassert>

#include "taichi/ir/ir_builder.h"
#include "taichi/ir/statements.h"
#include "taichi/program/program.h"
#include "taichi/ir/ir.h"

#include "taichi/gui/gui.h"

class TokenStream {
public:
    std::vector<std::string> tokens;
    std::vector<std::string>::iterator it;
    taichi::lang::IRBuilder *builder = nullptr;

    TokenStream(std::istream &in) {
        using namespace std;
        tokens.push_back("{");
        std::vector<std::string> &raw_tokens = tokens;
        for (;;) {
            for (; isspace(in.peek()); in.ignore());
            if (in.eof()) break;
            if (isalpha(in.peek())) {
                string tmp = "";
                for (;isalnum(in.peek()); tmp += in.get());
                raw_tokens.push_back(tmp);
            }
            else if (isdigit(in.peek())) {
                string tmp = "";
                for (;isdigit(in.peek()); tmp += in.get());
                raw_tokens.push_back(tmp);
            }
            else {
                string tmp = "";
                tmp += in.get();
                if (raw_tokens.size() && raw_tokens.back() == tmp) raw_tokens.back() += tmp;
                else raw_tokens.push_back(tmp);
            }
        }
        tokens.push_back("}");
        it = tokens.begin();
    }

    std::string& peek() {
        return *it;
    }

    std::string& get() {
        it++;
        return *(it - 1);
    }

    bool eof() {
        return it == tokens.end();
    }

    std::string& reset() {
        return *(it = tokens.begin());
    }
};

class CellNode {
public:
    std::string token;
    std::vector<CellNode*> children;
    static TokenStream *ts;

    CellNode(std::string _token): token(_token) {}

    static std::string get() {
        return ts->get();
    }

    static std::string peek() {
        return ts->peek();
    }

    static std::string eat(std::string str) {
        assert(get() == str);
        return str;
    }

    static CellNode* parse_single_element() {
        if (peek() == "(") {
            eat("(");
            CellNode *ans = parse_exp();
            eat(")");
            return ans;
        }
        if (isdigit(peek()[0])) return new CellNode(get());
        if (isalpha(peek()[0])) {
            CellNode *p = new CellNode(get());
            if (peek() == "[") {
                eat("[");
                p->children.push_back(parse_exp());
                for (; peek() == ",";) {
                    eat(",");
                    p->children.push_back(parse_exp());
                }
                eat("]");
            }
            return p;
        }
        CellNode *p = new CellNode(get());
        p->children.push_back(parse_single_element());
        return p;
    }

    static CellNode* parse_exp() {
        using namespace std;
        stack<CellNode*> eles;
        stack<string> ops;
        map<string, int> precedence;
        precedence["+"] = 6;
        precedence["-"] = 6;
        precedence["&"] = 11;
        precedence["|"] = 13;
        precedence["=="] = 10;
        precedence["||"] = 15;
        auto fun = [&]() {
            CellNode *first, *second;
            second = eles.top();
            eles.pop();
            first = eles.top();
            eles.pop();
            CellNode *p = new CellNode(ops.top());
            ops.pop();
            p->children.push_back(first);
            p->children.push_back(second);
            eles.push(p);
        };
        eles.push(parse_single_element());
        for (;peek() != ";" && peek() != ")" && peek() != "," && peek() != "]";) {
            assert(precedence.find(peek()) != precedence.end());
            for (; !ops.empty() && precedence[ops.top()] <= precedence[peek()]; fun());
            ops.push(get());
            eles.push(parse_single_element());
        }
        for (; !ops.empty(); fun());
        return eles.top();
    }

    static CellNode* parse_if() {
        CellNode *p = new CellNode(eat("if"));
        eat("(");
        p->children.push_back(parse_exp());
        eat(")");
        p->children.push_back(parse_stmt());
        if (peek() == "else") {
            eat("else");
            p->children.push_back(parse_stmt());
        }
        return p;
    }

    static CellNode* parse_stmt() {
        if (peek() == "{") return parse_block();
        if (peek() == "if") return parse_if();
        if (peek() == "int") {
            CellNode *p = new CellNode(eat("int"));
            p->children.push_back(new CellNode(get()));
            eat(";");
            return p;
        }
        else {
            CellNode *p = new CellNode("=");
            p->children.push_back(new CellNode(get()));
            eat("=");
            p->children.push_back(parse_exp());
            eat(";");
            return p;
        }
    }

    static CellNode* parse_block() {
        CellNode *p = new CellNode("{");
        eat("{");
        for (;peek() != "}"; p->children.push_back(parse_stmt()));
        eat("}");
        return p;
    }

    static CellNode* parse(TokenStream &_ts) {
        ts = &_ts;
        return parse_block();
    }
};