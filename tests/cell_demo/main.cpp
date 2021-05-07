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

void test_snode() {
    using namespace taichi;
    using namespace lang;
    auto prog_ = Program(arch_from_name("x64"));
    //prog_.thread_pool = std::make_unique<ThreadPool>(1);
    CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    prog_.config = config_print_ir;  // print ir
    auto root = prog_.snode_root.get();
    auto &dense = root->dense(Index(0), 10);
    auto &place = dense.insert_children(SNodeType::place);
    place.dt = PrimitiveType::i32;
    prog_.materialize_layout();


    IRBuilder builder;
    auto *zero = builder.get_int32(0);
    auto *ten = builder.get_int32(10);
    auto *loop = builder.create_range_for(zero, ten, 1, 0, 4);
    {
        builder.set_insertion_point_to_loop_begin(loop);
        auto *index = builder.get_loop_index(loop, 0);
        auto *ptr = builder.insert(Stmt::make<GlobalPtrStmt>(LaneAttribute<SNode *>(&place), std::vector<Stmt*>(1, index)));
        builder.insert(std::make_unique<GlobalStoreStmt>(ptr, index));
        builder.set_insertion_point_to_after(loop);
    }
    auto *index = builder.get_int32(4);
    auto *ret_ptr = builder.insert(Stmt::make<GlobalPtrStmt>(LaneAttribute<SNode *>(&place), std::vector<Stmt*>(1, index)));
    auto *ret_val = builder.insert(std::make_unique<GlobalLoadStmt>(ret_ptr));
    builder.create_return(ret_val);


    auto block = builder.extract_ir();
    auto ker = std::make_unique<Kernel>(prog_, std::move(block));
    auto launch_ctx = ker->make_launch_context();
    (*ker)(launch_ctx);
    std::cout << prog_.fetch_result<int>(0) << " : ans\n";
    exit(0);
}

int main() {
    test_snode();
    return 0;
}