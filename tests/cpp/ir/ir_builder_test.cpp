#include "gtest/gtest.h"

#include "taichi/ir/ir_builder.h"
#include "taichi/ir/statements.h"
#include "taichi/program/program.h"

namespace taichi {
namespace lang {

TEST(IRBuilder, Basic) {
  IRBuilder builder;
  auto *lhs = builder.get_int32(40);
  auto *rhs = builder.get_int32(2);
  auto *add = builder.create_add(lhs, rhs);
  ASSERT_TRUE(add->is<BinaryOpStmt>());
  auto *addc = add->cast<BinaryOpStmt>();
  EXPECT_EQ(addc->lhs, lhs);
  EXPECT_EQ(addc->rhs, rhs);
  EXPECT_EQ(addc->op_type, BinaryOpType::add);
  auto ir = builder.extract_ir();
  ASSERT_TRUE(ir->is<Block>());
  EXPECT_EQ(ir->as<Block>()->size(), 3);
}

TEST(IRBuilder, Print) {
  IRBuilder builder;
  auto *one = builder.get_int32(1);
  ASSERT_TRUE(one->is<ConstStmt>());
  std::string message = "message";
  auto *result = builder.create_print(one, message, one);
  ASSERT_TRUE(result->is<PrintStmt>());
  auto *print = result->cast<PrintStmt>();
  EXPECT_EQ(print->contents.size(), 3);
  ASSERT_TRUE(std::holds_alternative<Stmt *>(print->contents[0]));
  EXPECT_EQ(std::get<Stmt *>(print->contents[0]), one);
  ASSERT_TRUE(std::holds_alternative<std::string>(print->contents[1]));
  EXPECT_EQ(std::get<std::string>(print->contents[1]), message);
  ASSERT_TRUE(std::holds_alternative<Stmt *>(print->contents[2]));
  EXPECT_EQ(std::get<Stmt *>(print->contents[2]), one);
}

TEST(IRBuilder, RangeFor) {
  IRBuilder builder;
  auto *zero = builder.get_int32(0);
  auto *ten = builder.get_int32(10);
  auto *loop = builder.create_range_for(zero, ten);
  builder.set_insertion_point_to_loop_begin(loop);
  auto *index = builder.get_loop_index(loop, 0);
  builder.set_insertion_point_to_after(loop);
  auto *ret = builder.create_return(zero);
  EXPECT_EQ(zero->parent->size(), 4);
  ASSERT_TRUE(loop->is<RangeForStmt>());
  auto *loopc = loop->cast<RangeForStmt>();
  EXPECT_EQ(loopc->body->size(), 1);
  EXPECT_EQ(loopc->body->statements[0].get(), index);
}

TEST(IRBuilder, RunKernel) {
    auto prog_ = Program(arch_from_name("x64"));
    /*CompileConfig config_print_ir;
    config_print_ir.print_ir = true;
    prog_.config = config_print_ir;*/  // print ir
    prog_.materialize_layout();
    IRBuilder builder;
    auto *one = builder.get_int32(1);
    builder.create_return(one);
    auto block = builder.extract_ir();
    auto ker = std::make_unique<Kernel>(prog_, std::move(block));
    auto launch_ctx = ker->make_launch_context();
    (*ker)(launch_ctx);
    EXPECT_EQ(prog_.fetch_result<int64_t>(0), 1);
}

TEST(IRBuilder, RunSnodeKernel) {
    auto prog_ = Program(arch_from_name("x64"));
    prog_.materialize_layout();
    auto root = prog_.snode_root.get();
    auto &dense = root->dense(Index(0), 10);
    auto &place = dense.insert_children(SNodeType::place);
    place.dt = PrimitiveType::f32;
    IRBuilder builder;
    auto *index = builder.get_int32(1);
    //auto *ptr = builder.insert(std::make_unique<GlobalPtrStmt>(&dense, std::vector<Stmt*>(1, index)));
    auto *ptr = builder.insert(Stmt::make<GlobalPtrStmt>(LaneAttribute<SNode *>(&place), std::vector<Stmt*>(1, index)));
    //builder.insert(std::make_unique<GlobalStoreStmt>(ptr, index));
    builder.create_return(index);
    auto block = builder.extract_ir();
    auto ker = std::make_unique<Kernel>(prog_, std::move(block));
    auto launch_ctx = ker->make_launch_context();
    (*ker)(launch_ctx);
    EXPECT_EQ(prog_.fetch_result<int64_t>(0), 1);
}

}  // namespace lang
}  // namespace taichi
