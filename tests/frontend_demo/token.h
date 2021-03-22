#include <string>

class Token {
public:
    Token() {}
    virtual ~Token() = default;
    virtual std::string to_string() = 0;
};

class Constant : public Token {
public:
    double value;
    Constant(double v): Token(), value(v) {}
    std::string to_string() { return std::to_string(value); }
};

class Operator : public Token {
public:
    Operator(): Token() {}
    virtual int get_precedence() = 0;
};

class Add : public Operator {
public:
    Add(): Operator() {}
    int get_precedence() { return 13; }
    std::string to_string() { return "+"; }
};

class Sub : public Operator {
public:
    Sub(): Operator() {}
    int get_precedence() { return 13; }
    std::string to_string() { return "-"; }
};

class Mult : public Operator {
public:
    Mult(): Operator() {}
    int get_precedence() { return 14; }
    std::string to_string() { return "*"; }
};

class Div : public Operator {
public:
    Div(): Operator() {}
    int get_precedence() { return 14; }
    std::string to_string() { return "/"; }
};

class OpenParen : public Operator {
public:
    OpenParen(): Operator() {}
    int get_precedence() { return -1; }
    std::string to_string() { return "("; }
};

class CloseParen : public Operator {
public:
    CloseParen(): Operator() {}
    int get_precedence() { return -2; }
    std::string to_string() { return ")"; }
};

class Variable : public Token {
public:
    std::string name;
    Variable(std::string _name): Token(), name(_name) {}
    std::string to_string() { return name; }
};

class UnaryToken : public Operator {
public:
    std::string name;
    UnaryToken(std::string _name): Operator(), name(_name) {}
    int get_precedence() { return 15; }
    std::string to_string() { return isalpha(name[0]) ? "ti." + name : name; }
};
