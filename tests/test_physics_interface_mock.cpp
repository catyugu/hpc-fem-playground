/**
 * @file test_physics_interface_mock.cpp
 * @brief Mock test for PhysicsInterface contract verification
 * 

 * 
 * This test verifies the contract of the PhysicsInterface:
 * - Mock physics can be instantiated
 * - assemble() method can be called
 * - getFiniteElementSpace() returns valid pointer
 * - Interface polymorphism works correctly
 */

#include "hpcfem/core/physics_interface.hpp"
#include <gtest/gtest.h>
#include <memory>

using namespace hpcfem;

/**
 * @class MockPhysics
 * @brief Mock implementation of PhysicsInterface for testing
 */
class MockPhysics : public PhysicsInterface
{
public:
#ifdef MFEM_USE_MPI
    explicit MockPhysics(mfem::ParMesh* pmesh)
        : assembleCallCount(0)
    {
        fec = new mfem::H1_FECollection(1, pmesh->Dimension());
        fespace = new mfem::ParFiniteElementSpace(pmesh, fec);
    }
#else
    explicit MockPhysics(mfem::Mesh* mesh)
        : assembleCallCount(0)
    {
        fec = new mfem::H1_FECollection(1, mesh->Dimension());
        fespace = new mfem::FiniteElementSpace(mesh, fec);
    }
#endif

    ~MockPhysics() override
    {
        delete fespace;
        delete fec;
    }

#ifdef MFEM_USE_MPI
    void assemble(mfem::HypreParMatrix& A,
                 mfem::Vector& b,
                 mfem::Vector& x,
                 mfem::Array<int>& essTdofList) override
    {
        assembleCallCount++;
        // Mock implementation: create identity matrix
        int size = fespace->GetTrueVSize();
        b.SetSize(size);
        b = 1.0;
        x.SetSize(size);
        x = 0.0;
    }

    mfem::ParFiniteElementSpace* getFiniteElementSpace() override
    {
        return fespace;
    }
#else
    void assemble(mfem::SparseMatrix& A,
                 mfem::Vector& b,
                 mfem::Vector& x,
                 mfem::Array<int>& essTdofList) override
    {
        assembleCallCount++;
        // Mock implementation: create identity matrix
        int size = fespace->GetTrueVSize();
        b.SetSize(size);
        b = 1.0;
        x.SetSize(size);
        x = 0.0;
    }

    mfem::FiniteElementSpace* getFiniteElementSpace() override
    {
        return fespace;
    }
#endif

    int getAssembleCallCount() const { return assembleCallCount; }

private:
    mfem::H1_FECollection* fec;
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* fespace;
#else
    mfem::FiniteElementSpace* fespace;
#endif
    int assembleCallCount;
};

/**
 * @brief Test that MockPhysics implements PhysicsInterface correctly
 */
TEST(PhysicsInterfaceTest, MockPhysicsInstantiation)
{
    constexpr int MESH_SIZE = 2;
    
#ifdef MFEM_USE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian3D(
        MESH_SIZE, MESH_SIZE, MESH_SIZE,
        mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0
    );
    mfem::ParMesh pmesh(comm, serialMesh);
    MockPhysics physics(&pmesh);
#else
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
        MESH_SIZE, MESH_SIZE, MESH_SIZE,
        mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0
    );
    MockPhysics physics(&mesh);
#endif
    
    EXPECT_EQ(physics.getAssembleCallCount(), 0);
    EXPECT_NE(physics.getFiniteElementSpace(), nullptr);
}

/**
 * @brief Test that assemble() method can be called through interface
 */
TEST(PhysicsInterfaceTest, AssembleMethodCall)
{
    constexpr int MESH_SIZE = 2;
    
#ifdef MFEM_USE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian3D(
        MESH_SIZE, MESH_SIZE, MESH_SIZE,
        mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0
    );
    mfem::ParMesh pmesh(comm, serialMesh);
    MockPhysics physics(&pmesh);
    
    mfem::HypreParMatrix A;
#else
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
        MESH_SIZE, MESH_SIZE, MESH_SIZE,
        mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0
    );
    MockPhysics physics(&mesh);
    
    mfem::SparseMatrix A;
#endif
    
    mfem::Vector b, x;
    mfem::Array<int> essTdofList;
    
    // Call assemble through interface
    PhysicsInterface* physicsPtr = &physics;
    physicsPtr->assemble(A, b, x, essTdofList);
    
    // Verify assemble was called
    EXPECT_EQ(physics.getAssembleCallCount(), 1);
    
    // Verify vectors were sized
    EXPECT_GT(b.Size(), 0);
    EXPECT_GT(x.Size(), 0);
}

/**
 * @brief Test polymorphism through interface pointer
 */
TEST(PhysicsInterfaceTest, Polymorphism)
{
    constexpr int MESH_SIZE = 2;
    
#ifdef MFEM_USE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian3D(
        MESH_SIZE, MESH_SIZE, MESH_SIZE,
        mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0
    );
    mfem::ParMesh pmesh(comm, serialMesh);
    std::unique_ptr<PhysicsInterface> physics = std::make_unique<MockPhysics>(&pmesh);
    
    mfem::HypreParMatrix A;
#else
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
        MESH_SIZE, MESH_SIZE, MESH_SIZE,
        mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0
    );
    std::unique_ptr<PhysicsInterface> physics = std::make_unique<MockPhysics>(&mesh);
    
    mfem::SparseMatrix A;
#endif
    
    EXPECT_NE(physics.get(), nullptr);
    EXPECT_NE(physics->getFiniteElementSpace(), nullptr);
    
    mfem::Vector b, x;
    mfem::Array<int> essTdofList;
    
    // Verify we can call assemble through base class pointer
    physics->assemble(A, b, x, essTdofList);
    
    EXPECT_GT(b.Size(), 0);
}

/**
 * @brief Test that getFiniteElementSpace returns valid space
 */
TEST(PhysicsInterfaceTest, FiniteElementSpace)
{
    constexpr int MESH_SIZE = 2;
    
#ifdef MFEM_USE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian3D(
        MESH_SIZE, MESH_SIZE, MESH_SIZE,
        mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0
    );
    mfem::ParMesh pmesh(comm, serialMesh);
    MockPhysics physics(&pmesh);
    
    mfem::ParFiniteElementSpace* fespace = physics.getFiniteElementSpace();
#else
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
        MESH_SIZE, MESH_SIZE, MESH_SIZE,
        mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0
    );
    MockPhysics physics(&mesh);
    
    mfem::FiniteElementSpace* fespace = physics.getFiniteElementSpace();
#endif
    
    EXPECT_NE(fespace, nullptr);
    EXPECT_GT(fespace->GetTrueVSize(), 0);
    EXPECT_EQ(fespace->GetVDim(), 1);
}

int main(int argc, char** argv)
{
#ifdef MFEM_USE_MPI
    MPI_Init(&argc, &argv);
#endif
    
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    
#ifdef MFEM_USE_MPI
    MPI_Finalize();
#endif
    
    return result;
}
