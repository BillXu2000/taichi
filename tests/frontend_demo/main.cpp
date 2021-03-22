#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <stack>
#include <set>
#include <cassert>
#include "expr.h"

std::vector<std::unique_ptr<Token> > get_lexical(std::string plain, std::set<std::string> &names, std::set<std::string> reserved_funcs) {
    using namespace std;
    vector<unique_ptr<Token> > lex;
    auto is_number_or_dot = [](char ch) {
        return isdigit(ch) || ch == '.';
    };
    bool expect_operator = false;
    lex.push_back(make_unique<OpenParen>());
    for (int i = 0; plain[i];) {
        if (plain[i] == ';') break;
        if (isspace(plain[i])) i++;
        else if (is_number_or_dot(plain[i])) {
            assert(!expect_operator);
            expect_operator = true;
            int j;
            for (j = i + 1; is_number_or_dot(plain[j]); j++);
            lex.push_back(make_unique<Constant>(atof(plain.substr(i, j - i).c_str())));
            i = j;
        }
        else if (isalpha(plain[i])) {
            assert(!expect_operator);
            int j;
            for (j = i + 1; isalpha(plain[j]); j++);
            string name = plain.substr(i, j - i);
            if (reserved_funcs.find(name) == reserved_funcs.end()) {
                lex.push_back(make_unique<Variable>(name));
                names.insert(name);
                expect_operator = true;
            }
            else lex.push_back(make_unique<UnaryToken>(name));
            i = j;
        }
        else {
            if (plain[i] == '-') {
                if (expect_operator) lex.push_back(make_unique<Sub>());
                else lex.push_back(make_unique<UnaryToken>("-"));
                expect_operator = false;
            }
            else if (plain[i] == '(') {
                assert(!expect_operator);
                lex.push_back(make_unique<OpenParen>());
            }
            else if (plain[i] == ')') {
                assert(expect_operator);
                lex.push_back(make_unique<CloseParen>());
            }
            else {
                assert(expect_operator);
                expect_operator = false;
                switch (plain[i]) {
                    case '+':
                        lex.push_back(make_unique<Add>());
                        break;
                    case '*':
                        lex.push_back(make_unique<Mult>());
                        break;
                    case '/':
                        lex.push_back(make_unique<Div>());
                        break;
                }
            }
            i++;
        }
    }
    lex.push_back(make_unique<CloseParen>());
    return lex;
}

Expr* get_ast(std::vector<std::unique_ptr<Token> >& lex) {
    using namespace std;
    stack<unique_ptr<Operator> > ops;
    stack<unique_ptr<Expr> > nodes;
    for (auto &i : lex) {
        auto op = dynamic_cast<Operator*>(i.get());
        if (op) {
            if (dynamic_cast<OpenParen*>(op) || dynamic_cast<UnaryToken*>(op)) {
                i.release();
                ops.push(unique_ptr<Operator>(op));
                continue;
            }
            for (;!ops.empty() && ops.top()->get_precedence() >= op->get_precedence();) {
                if (dynamic_cast<OpenParen*>(ops.top().get())) {
                    ops.pop();
                    break;
                }
                auto right = nodes.top().get();
                nodes.top().release();
                nodes.pop();
                if (dynamic_cast<UnaryToken*>(ops.top().get())) {
                    nodes.push(make_unique<UnaryOp>(move(ops.top()), unique_ptr<Expr>(right)));
                }
                else {
                    auto left = nodes.top().get();
                    nodes.top().release();
                    nodes.pop();
                    nodes.push(make_unique<BinOp>(unique_ptr<Expr>(left), move(ops.top()), unique_ptr<Expr>(right)));
                }
                ops.pop();
            }
            if (!dynamic_cast<CloseParen*>(op)) {
                i.release();
                ops.push(unique_ptr<Operator>(op));
            }
        }
        else nodes.push(make_unique<Literal>(move(i)));
    }
    auto ans = nodes.top().get();
    nodes.top().release();
    return ans;
}

void print(std::unique_ptr<Expr> &i) {
    using namespace std;
    auto p = i.get();
    auto bin = dynamic_cast<BinOp*>(p);
    auto unary = dynamic_cast<UnaryOp*>(p);
    if (bin) {
        cout << bin->op->to_string() << "(";
        print(bin->left);
        cout << ", ";
        print(bin->right);
        cout << ")";
    }
    else if (unary) {
        cout << unary->op->to_string() << "(";
        print(unary->operand);
        cout << ")";
    }
    else cout << dynamic_cast<Literal*>(p)->value->to_string();
}

void init_IR(std::fstream &fs, std::set<std::string> &names) {
    fs << ""
"import taichi as ti\n"
//"n_IR = 100000\n"
"real = ti.f32\n";
//"ssa = ti.field(dtype = real, shape = n_IR, needs_grad = True)\n";
//"ans = ti.field(dtype = real, shape = (), needs_grad = True)\n";
    for (auto i : names) {
        if (i != "x" && i != "y") fs << i + " = ti.field(dtype = real, shape = (), needs_grad = True)\n";
    }
    fs << ""
"@ti.func\n"
"def generated_ir(x, y):\n";
}

void build_IR(std::unique_ptr<Expr> &i, std::fstream &fs, int &cnt) {
    using namespace std;
    auto p = i.get();
    auto bin = dynamic_cast<BinOp*>(p);
    auto unary = dynamic_cast<UnaryOp*>(p);
    if (bin) {
        build_IR(bin->left, fs, cnt);
        build_IR(bin->right, fs, cnt);
        p->ssa = cnt++;
        //fs << "\tssa[" + to_string(bin->ssa) + "] = ssa[" + to_string(bin->left->ssa) + "] " + bin->op->to_string() + " ssa[" + to_string(bin->right->ssa) + "]\n";
        fs << "\tssa" + to_string(bin->ssa) + " = ssa" + to_string(bin->left->ssa) + " " + bin->op->to_string() + " ssa" + to_string(bin->right->ssa) + "\n";
    }
    else if (unary) {
        build_IR(unary->operand, fs, cnt);
        p->ssa = cnt++;
        //fs << "\tssa[" + to_string(unary->ssa) + "] = " + unary->op->to_string() + "(ssa[" + to_string(unary->operand->ssa) + "])\n";
        fs << "\tssa" + to_string(unary->ssa) + " = " + unary->op->to_string() + "(ssa" + to_string(unary->operand->ssa) + ")\n";
    }
    else {
        p->ssa = cnt++;
        //fs << "\tssa[" + to_string(p->ssa) + "] = " + dynamic_cast<Literal*>(p)->value->to_string() + "\n";
        fs << "\tssa" + to_string(p->ssa) + " = " + dynamic_cast<Literal*>(p)->value->to_string() + "\n";
    }
}


int main(int argc, char* argv[]) {
    std::fstream fs(argv[1], std::fstream::in);
    std::string plain;
    for (char ch; fs.get(ch); plain += ch);
    fs.close();
    std::set<std::string> names, reserved_funcs;
    reserved_funcs.insert("cos");
    reserved_funcs.insert("sin");
    auto &&lex = get_lexical(plain, names, reserved_funcs);
    //printf("plain: %d, lex: %d\n", plain.size(), lex.size());
    puts("before");
    auto root = std::unique_ptr<Expr>(get_ast(lex));
    puts("after");
    print(root);
    puts("");
    fs.open(argv[2], std::fstream::out);
    init_IR(fs, names);
    int cnt = 0;
    build_IR(root, fs, cnt);
    //fs << "\treturn ssa[" + std::to_string(cnt - 1) + "]";
    fs << "\treturn ssa" + std::to_string(cnt - 1) + "";
    return 0;
}
