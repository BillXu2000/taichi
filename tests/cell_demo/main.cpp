#include "parser.h"

class BuilderHelper;

typedef BuilderHelper BH;

class BuilderHelper {
public:
    static taichi::lang::IRBuilder *builder;
    static std::map<std::string, taichi::lang::AllocaStmt*> symbol;
    static std::map<std::string, taichi::lang::SNode*> snode_table;
    static std::vector<taichi::lang::Stmt*> global_indices;
    taichi::lang::Stmt *stmt;

    BuilderHelper(taichi::lang::Stmt *_stmt): stmt(_stmt) {}

    BuilderHelper operator + (taichi::lang::Stmt *other) {
        return builder->create_add(stmt, other);
    }

    BuilderHelper operator - (taichi::lang::Stmt *other) {
        return builder->create_sub(stmt, other);
    }

    BuilderHelper operator * (taichi::lang::Stmt *other) {
        return builder->create_mul(stmt, other);
    }

    BuilderHelper operator / (taichi::lang::Stmt *other) {
        return builder->create_div(stmt, other);
    }

    BuilderHelper operator % (taichi::lang::Stmt *other) {
        return builder->create_mod(stmt, other);
    }

    BuilderHelper operator > (taichi::lang::Stmt *other) {
        return builder->create_cmp_gt(stmt, other);
    }

    BuilderHelper operator < (taichi::lang::Stmt *other) {
        return builder->create_cmp_lt(stmt, other);
    }

    BuilderHelper operator | (taichi::lang::Stmt *other) {
        return builder->create_or(stmt, other);
    }

    BuilderHelper operator == (taichi::lang::Stmt *other) {
        return builder->create_cmp_eq(stmt, other);
    }

    BuilderHelper operator & (taichi::lang::Stmt *other) {
        return builder->create_and(stmt, other);
    }

    operator taichi::lang::Stmt*() {
        return stmt;
    }

    static BH gen(CellNode *p) {
        using namespace std;
        using namespace taichi;
        using namespace lang;
        if (p->token == "{") {
            for (auto i : p->children) gen(i);
            return nullptr;
        }
        if (p->token == "int") {
            symbol[p->children[0]->token] = builder->create_local_var(PrimitiveType::i32);
            return nullptr;
        }
        if (p->token == "=") {
            builder->create_local_store(symbol[p->children[0]->token], gen(p->children[1]));
            return nullptr;
        }
        if (p->token == "if") {
            auto *branch = builder->create_if(gen(p->children[0]));
            {
                auto _ = builder->get_if_guard(branch, true);
                gen(p->children[1]);
            }
            if (p->children.size() > 2) {
                auto _ = builder->get_if_guard(branch, false);
                gen(p->children[2]);
            }
            return nullptr;
        }
        if (isdigit(p->token[0])) return builder->get_int32(atoi(p->token.c_str()));
        if (isalpha(p->token[0])) {
            if (p->children.size()) {
                vector<Stmt*> indices;
                for (int i = 0; i < p->children.size(); i++) {
                    indices.push_back(BH(global_indices[i]) + gen(p->children[i]));
                }
                return builder->create_global_load(builder->create_global_ptr(snode_table[p->token], indices));
            }
            return builder->create_local_load(symbol[p->token]);
        }
        if (p->children.size() == 1) {
            if (p->token == "-") return builder->create_neg(gen(p->children[0]));
        }
        if (p->children.size() == 2) {
            if (p->token == "+") return builder->create_add(gen(p->children[0]), gen(p->children[1]));
            if (p->token == "-") return builder->create_sub(gen(p->children[0]), gen(p->children[1]));
            if (p->token == "&") return builder->create_and(gen(p->children[0]), gen(p->children[1]));
            if (p->token == "|") return builder->create_or(gen(p->children[0]), gen(p->children[1])); // same as ||
            if (p->token == "||") return builder->create_or(gen(p->children[0]), gen(p->children[1]));
            if (p->token == "==") return builder->create_cmp_eq(gen(p->children[0]), gen(p->children[1]));
        }
        assert(false);
    }
};
taichi::lang::IRBuilder *BuilderHelper::builder(nullptr);
std::map<std::string, taichi::lang::AllocaStmt*> BH::symbol;
std::map<std::string, taichi::lang::SNode*> BH::snode_table;
std::vector<taichi::lang::Stmt*> BH::global_indices;

class BuilderHelperGuard {
public:
    BuilderHelperGuard(taichi::lang::IRBuilder &builder) {
        BuilderHelper::builder = &builder;
    }

    ~BuilderHelperGuard() {
        BuilderHelper::builder = nullptr;
    }
};

void game_of_life(CellNode *cell_root, int argc, char *argv[]) {
    using namespace std;
    using namespace taichi;
    using namespace lang;
    auto program = Program(arch_from_name("x64"));
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    program.config = config_print_ir;  // print ir
    auto root = program.snode_root.get();

    const int N = 512, M = 16;
    std::vector<Index> index_dense = {0, 1};
    std::vector<int> size_dense = {M, M};
    auto &pointer = root->pointer(vector<Index>{0, 1}, vector<int>{N / M, N / M}).dense(index_dense, size_dense);
    auto &alive = pointer.insert_children(SNodeType::place);
    //auto &alive = root->dense(index_dense, size_dense).insert_children(SNodeType::place);
    alive.dt = PrimitiveType::i32;
    auto &next = root->pointer(vector<Index>{0, 1}, vector<int>{N / M, N / M}).dense(index_dense, size_dense).insert_children(SNodeType::place);
    //auto &next = root->dense(index_dense, size_dense).insert_children(SNodeType::place);
    next.dt = PrimitiveType::i32;

    program.materialize_layout();

    typedef BuilderHelper BH;

    /*std::unique_ptr<Kernel> kernel_init;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *left = builder.get_int32(M / 2 * N);
        auto *right = builder.get_int32(M / 2 * N + N);
        auto *loop_x = builder.create_range_for(left, right, 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop_x);
            auto *x = builder.get_loop_index(loop_x, 0);
            auto *loop_y = builder.create_range_for(left, right, 1, 0, 4);
            {
                auto _ = builder.get_loop_guard(loop_y);
                auto *y = builder.get_loop_index(loop_y, 0);
                Stmt *cnt = BH(x) * builder.get_int32(N) + y;
                auto *mod = builder.get_int32(10007);
                //Stmt *ans = (BH(cnt) * builder.get_int32(9287) % mod + mod) % builder.get_int32(2);
                //Stmt *rnd = ((BH(cnt) * builder.get_int32(9287) % mod + mod) % builder.get_int32(16));
                //Stmt *rnd = BH(builder.insert(make_unique<RandStmt>(PrimitiveType::i32))) % builder.get_int32(16);
                Stmt *rnd = BH(builder.create_rand(PrimitiveType::i32)) % builder.get_int32(2);
                Stmt *ans = (BH(x) > builder.get_int32(M / 2 * N + N / 2) | BH(y) > builder.get_int32(M / 2 * N + N / 2)) & rnd & (BH(rnd) < builder.get_int32(16));
                vector<Stmt *> indices = {x, y};
                builder.create_global_store(builder.create_global_ptr(&alive, indices), ans);
            }
        }
        kernel_init = make_unique<Kernel>(program, builder.extract_ir(), "init");
    }*/

    std::unique_ptr<Kernel> kernel_init;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        if (argc > 1) {
            fstream in(argv[1], fstream::in);
            int n, m;
            in >> n >> m;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    int k;
                    in >> k;
                    vector<Stmt*> indices{builder.get_int32(i + 2), builder.get_int32(j + 2)};
                    if (k) builder.create_global_store(builder.create_global_ptr(&alive, indices), builder.get_int32(1));
                }
            }
        }
        else {
            auto *left = builder.get_int32(M / 2 * N);
            auto *right = builder.get_int32(M / 2 * N + N);
            auto *loop_x = builder.create_range_for(left, right, 1, 0, 4);
            {
                auto _ = builder.get_loop_guard(loop_x);
                auto *x = builder.get_loop_index(loop_x, 0);
                auto *loop_y = builder.create_range_for(left, right, 1, 0, 4);
                {
                    auto _ = builder.get_loop_guard(loop_y);
                    auto *y = builder.get_loop_index(loop_y, 0);
                    Stmt *cnt = BH(x) * builder.get_int32(N) + y;
                    auto *mod = builder.get_int32(10007);
                    //Stmt *ans = (BH(cnt) * builder.get_int32(9287) % mod + mod) % builder.get_int32(2);
                    //Stmt *rnd = ((BH(cnt) * builder.get_int32(9287) % mod + mod) % builder.get_int32(16));
                    //Stmt *rnd = BH(builder.insert(make_unique<RandStmt>(PrimitiveType::i32))) % builder.get_int32(16);
                    Stmt *rnd = BH(builder.create_rand(PrimitiveType::i32)) % builder.get_int32(2);
                    Stmt *ans = (BH(x) > builder.get_int32(M / 2 * N + N / 2) | BH(y) > builder.get_int32(M / 2 * N + N / 2)) & rnd & (BH(rnd) < builder.get_int32(16));
                    vector<Stmt *> indices = {x, y};
                    builder.create_global_store(builder.create_global_ptr(&alive, indices), ans);
                }
            }
        }
        kernel_init = make_unique<Kernel>(program, builder.extract_ir(), "init");
    }
    
    std::unique_ptr<Kernel> kernel_step;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *loop = builder.create_struct_for(&pointer, -1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop);
            BH x = builder.get_loop_index(loop, 0);
            BH y = builder.get_loop_index(loop, 1);
            BH::snode_table["alive"] = &alive;
            BH::global_indices.push_back(x);
            BH::global_indices.push_back(y);
            BH::symbol["output"] = builder.create_local_var(PrimitiveType::i32);
            BH::gen(cell_root);
            vector<Stmt *> indices = {x, y};
            //builder.create_global_store(builder.create_global_ptr(&next, indices), builder.create_local_load(BH::symbol["output"]));
            BH self = builder.create_global_load(builder.create_global_ptr(&alive, indices));
            BH ans = builder.create_local_load(BH::symbol["output"]);
            /*ans = ans + ((BH(x) == builder.get_int32(M / 2 * N + N - 1)) & builder.get_int32(1) & ((BH(self) & builder.get_int32(4)) > builder.get_int32(0)));
            ans = ans + ((BH(y) == builder.get_int32(M / 2 * N + N - 1)) & builder.get_int32(2) & ((BH(self) & builder.get_int32(8)) > builder.get_int32(0)));
            ans = ans + ((BH(x) == builder.get_int32(M / 2 * N)) & builder.get_int32(4) & ((BH(self) & builder.get_int32(1)) > builder.get_int32(0)));
            ans = ans + ((BH(y) == builder.get_int32(M / 2 * N)) & builder.get_int32(8) & ((BH(self) & builder.get_int32(2)) > builder.get_int32(0)));*/
            auto* branch = builder.create_if(self);
            {
                auto _ = builder.get_if_guard(branch, true);
                for (int i = -1; i < 2; i++) {
                    for (int j = -1; j < 2; j++) {
                        vector<Stmt *> indices = {x + builder.get_int32(i), y + builder.get_int32(j)};
                        //builder.insert(make_unique<AtomicOpStmt>(AtomicOpType::add, builder.create_global_ptr(&alive, indices), builder.get_int32(0)));
                        builder.create_atomic_add(builder.create_global_ptr(&alive, indices), builder.get_int32(0));
                    }
                }
            }
            builder.create_global_store(builder.create_global_ptr(&next, indices), ans);
            /*Stmt *sum = nullptr, *self = nullptr;
            for (int i = -1; i < 2; i++) {
                for (int j = -1; j < 2; j++) {
                    vector<Stmt *> indices = {
                        builder.create_add(builder.get_int32(i), x),
                        builder.create_add(builder.get_int32(j), y)};
                    auto *load = builder.create_global_load(builder.create_global_ptr(&alive, indices));
                    if (sum) sum = builder.create_add(sum, load);
                    else sum = load;
                    if (!i && !j) self = load;
                }
            }
            auto ans = (BH(sum) == builder.get_int32(3) | (BH(sum) - self) == builder.get_int32(3)) & builder.get_int32(1);
            vector<Stmt *> indices = {x, y};
            builder.create_global_store(builder.create_global_ptr(&next, indices), ans);*/
            /*BH ans = builder.get_int32(0);
            int way[4][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
            for (int i = 0; i < 4; i++) {
                vector<Stmt*> index = {way[i][0] ? BH(x) + builder.get_int32(way[i][0]) : x, way[i][1] ? BH(y) + builder.get_int32(way[i][1]) : y};
                BH load = builder.create_global_load(builder.create_global_ptr(&alive, index));
                ans = ans | load & builder.get_int32(1 << i);
            }
            vector<Stmt *> indices = {x, y};
            BH self = builder.create_global_load(builder.create_global_ptr(&alive, indices));
            ans = ans + ((BH(x) == builder.get_int32(N - 2)) & builder.get_int32(1) & ((BH(self) & builder.get_int32(4)) > builder.get_int32(0)));
            ans = ans + ((BH(y) == builder.get_int32(N - 2)) & builder.get_int32(2) & ((BH(self) & builder.get_int32(8)) > builder.get_int32(0)));
            ans = ans + ((BH(x) == builder.get_int32(1)) & builder.get_int32(4) & ((BH(self) & builder.get_int32(1)) > builder.get_int32(0)));
            ans = ans + ((BH(y) == builder.get_int32(1)) & builder.get_int32(8) & ((BH(self) & builder.get_int32(2)) > builder.get_int32(0)));
            ans = ans + ((ans == builder.get_int32(5)) & builder.get_int32(5)) - ((ans == builder.get_int32(10)) & builder.get_int32(5));
            builder.create_global_store(builder.create_global_ptr(&next, indices), ans);*/
        }
        kernel_step = make_unique<Kernel>(program, builder.extract_ir(), "step");
    }

    std::unique_ptr<Kernel> kernel_swap;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *loop = builder.create_struct_for(&pointer, 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop);
            auto *x = builder.get_loop_index(loop, 0);
            auto *y = builder.get_loop_index(loop, 1);
            vector<Stmt *> indices = {x, y};
            builder.create_global_store(
                builder.create_global_ptr(&alive, indices),
                builder.create_global_load(
                    builder.create_global_ptr(&next, indices)));
        }
        kernel_swap = make_unique<Kernel>(program, builder.extract_ir(), "swap");
    }

    /*std::unique_ptr<Kernel> kernel_gui;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *left = builder.get_int32(M / 2 * N);
        auto *right = builder.get_int32(M / 2 * N + N);
        auto *loop_x = builder.create_range_for(left, right, 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop_x);
            auto *x = builder.get_loop_index(loop_x, 0);
            auto *loop_y = builder.create_range_for(left, right, 1, 0, 4);
            {
                auto _ = builder.get_loop_guard(loop_y);
                auto *y = builder.get_loop_index(loop_y, 0);
                vector<Stmt *> indices = {x, y};
                builder.create_global_store(
                    builder.create_external_ptr(
                        builder.create_arg_load(0, PrimitiveType::i32, true),
                        vector<Stmt *>(1, BH(x) % builder.get_int32(N) * builder.get_int32(N) + BH(y) % builder.get_int32(N))),
                    builder.create_global_load(
                        builder.create_global_ptr(&alive, indices)));
            }
        }
        kernel_gui = make_unique<Kernel>(program, builder.extract_ir(), "gui");
        kernel_gui->insert_arg(PrimitiveType::gen, true);
    }*/
    std::unique_ptr<Kernel> kernel_gui;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *loop = builder.create_struct_for(&pointer, 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop);
            auto *x = builder.get_loop_index(loop, 0);
            auto *y = builder.get_loop_index(loop, 1);
            vector<Stmt *> indices = {x, y};
            builder.create_global_store(
                builder.create_external_ptr(
                    builder.create_arg_load(0, PrimitiveType::i32, true),
                    vector<Stmt *>(1, BH(x) % builder.get_int32(N) * builder.get_int32(N) + BH(y) % builder.get_int32(N))),
                BH(builder.create_global_load(builder.create_global_ptr(&alive, indices))) + builder.get_int32(2));
                //builder.create_global_load(builder.create_global_ptr(&alive, indices)));
        }
        kernel_gui = make_unique<Kernel>(program, builder.extract_ir(), "gui");
        kernel_gui->insert_arg(PrimitiveType::gen, true);
    }
    auto gui = GUI("cell language", N, N, true, false, 0, false, false);
    auto &canvas = *gui.canvas;
    static int ans[N][N];
    auto ctx_init = kernel_init->make_launch_context();
    auto ctx_step = kernel_step->make_launch_context();
    auto ctx_swap = kernel_swap->make_launch_context();
    auto ctx_gui = kernel_gui->make_launch_context();
    ctx_gui.set_arg_external_array(0, taichi::uint64(ans), 0);
    (*kernel_init)(ctx_init);
    for (int frame = 0;; frame++) {
        gui.update();
        //if (frame % 5) continue;
        (*kernel_step)(ctx_step);
        if (frame) (*kernel_swap)(ctx_swap);
        (*kernel_gui)(ctx_gui);
        long long sum = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                float k = ans[i][j] & 1, z = (ans[i][j] & 2) * 0.2;
                sum += ans[i][j];
                std::array<taichi::real, 4> color{k, z, z, 1};
                canvas.img[i][j] = taichi::Vector4(color);
                /*float k = 0;
                for (int z = 0; z < 4; z++) {
                    if ((ans[i][j] >> z) & 1) k += 0.25;
                }
                std::array<taichi::real, 4> color{k, 0, 0, 1};
                canvas.img[i][j] = taichi::Vector4(color);*/
            }
        }
        string name = to_string(frame);
        for (; name.size() < 5; name = "0" + name);
        if (frame % 59 == 0 || frame < 60) gui.screenshot(name + ".png");
    }
}

TokenStream *CellNode::ts = nullptr;

void cell_print(CellNode *p, int indent = 0) {
    for (int i = 0; i < indent; i++) {
        std::cout << " ";
    }
    std::cout << p->token << "\n";
    for (auto i : p->children) {
        cell_print(i, indent + 4);
    }
}

void test_pointer_struct_for() {
    using namespace std;
    using namespace taichi;
    using namespace lang;
    auto program = Program(arch_from_name("x64"));
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    program.config = config_print_ir;  // print ir
    auto root = program.snode_root.get();

    //auto &alive = root->pointer(vector<Index>{0}, vector<int>{9}).insert_children(SNodeType::place);
    auto &pointer = root->pointer(vector<Index>{0}, vector<int>{9});
    auto &alive = pointer.insert_children(SNodeType::place);
    alive.dt = PrimitiveType::i32;

    program.materialize_layout();

    typedef BuilderHelper BH;
    std::unique_ptr<Kernel> kernel_for;
    {
        IRBuilder builder;
        BuilderHelperGuard _(builder);
        auto *loop = builder.create_struct_for(&pointer, 1, 0, 4);
        //auto *loop = builder.create_range_for(builder.get_int32(0), builder.get_int32(9), 1, 0, 4);
        {
            auto _ = builder.get_loop_guard(loop);
            auto *x = builder.get_loop_index(loop, 0);
            builder.create_global_store(builder.create_global_ptr(&alive, vector<Stmt*>{x}), x);
        }
        kernel_for = make_unique<Kernel>(program, builder.extract_ir(), "pointer struct for");
    }
    auto ctx = kernel_for->make_launch_context();
    (*kernel_for)(ctx);
    exit(0);
}

int main(int argc, char *argv[]) {
    //test_pointer_struct_for();
    TokenStream ts(std::cin);
    CellNode *root = CellNode::parse(ts);
    cell_print(root);
    game_of_life(root, argc, argv);
    return 0;
}