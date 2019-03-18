#pragma once
#include <vector>
#include <bitset>
#include <queue>

#include "zjw_macro.h"
#include "util\zjw_math.h"
#include "util\zjw_obj.h"
#include "zjw_octree_double_buffer.h"


////////////////////////////////////////////////////////////////////////////////
//节点包含点的数据
////////////////////////////////////////////////////////////////////////////////

struct DBufferNodeData
{
	DBufferNodeData()
	{
	}

	~DBufferNodeData()
	{

	}
	//每个子空间有一堆的 particle的位置信息
	vector<Vec3> pointPosList;

	//每个在空间有一堆的 particle在obj中的序号
	vector<int> pointIdxList;

}; 

////////////////////////////////////////////////////////////////////////////////
//对target frame 或者是reference frame的相关内容进行删除
////////////////////////////////////////////////////////////////////////////////

class TraverseAllNodeForFree : public DoubleBufferOctree<DBufferNodeData>::Callback
{
public:
	bool isTarget;
public:
	TraverseAllNodeForFree(bool isTarget);
	~TraverseAllNodeForFree();

	/*
	*作用： 清理这个tree上的需要删除的node信息
	*返回参数: 返回false,终止递归遍历;返回true,继续递归遍历子节点
	*/
	virtual bool operator()(const Vec3 min, const Vec3 max, OctreeDoubelBufferNode<DBufferNodeData>* currNode);
};

////////////////////////////////////////////////////////////////////////////////
//遍历整个八叉树,处理叶子节点上的信息
////////////////////////////////////////////////////////////////////////////////

class TraverseGetInfoSetLeaf : public DoubleBufferOctree<DBufferNodeData>::Callback
{
public:
	//保存叶子节点的boundary ，用于渲染所用
	vector<Vec3> minVList;
	vector<Vec3> maxVList;

	//保存叶子节点的list,并给叶子节点编号
	vector<OctreeDoubelBufferNode<DBufferNodeData> *> nodeList;
	//叶子节点的编号器
	int leafIncr;
	//default : true
	bool isTarget;

public:
	TraverseGetInfoSetLeaf(bool isTarget_in = true);
	~TraverseGetInfoSetLeaf();

	void setIsTarget(bool flag);
	void initParam();

	virtual bool operator()(const Vec3 min, const Vec3 max, DBufferNodeData& nodeData) {
		return true;
	}
	//得到八叉树的叶子信息，把相应的信息保存到 参数中
	virtual bool operator()(const Vec3 min, const Vec3 max, OctreeDoubelBufferNode<DBufferNodeData>* currNode);
};

////////////////////////////////////////////////////////////////////////////////
//双buffer 八叉树
////////////////////////////////////////////////////////////////////////////////
// define frame identifier



class PcsDBufferOctree
{

public:
	//bounding box and cell size
	Vec3 minPos;
	Vec3 maxPos;
	Vec3 cellSize;

	//target(p-frame) 叶子节点的个数
	int refre_leaf_count_;
	int target_leaf_count_;

	//I frame 和P frame的指针，不用删除内容
	vector<Vertex> * referenceMeshVerList;
	vector<Vertex> * targetMeshVerList;

	//点云帧中的Point的数量
	int refFramePointNum;
	int targetFramePointNum;

	//当前八叉树中的具体情况
	bool targetFrameLoaded;
	bool referenceFrameLoaded;

	//双Buffer的序号标志 0 or 1
	int targetFrameId;
	int refrenceFrameId;

	DoubleBufferOctree<DBufferNodeData> *pcsDBufferOct;
	TraverseGetInfoSetLeaf* ctLeaf;
	//TraverseGetInfoSetLeaf *ctLeaf;
	
	
private:
	void setTargetFrameId(int id);
	void setReferenceFrameId(int id);
	
public:
	PcsDBufferOctree(Vec3 minPos,Vec3 maxPos, Vec3 cellSize);
	PcsDBufferOctree(Vec3 minPos, Vec3 maxPos, double cellSize);
	PcsDBufferOctree();

	~PcsDBufferOctree();

	/**拿到frame中的点的个数*/
	uint32_t getFramePointNum(bool isTarget = false);
	/**拿到octree的边界区*/
	void getBoundingBox(double &min_x, double &min_y, double &min_z, double & max_x, double & max_y, double &max_z);
	void setBoundingBox(const double &min_x, const double &min_y, const double &min_z, 
		const double & max_x, const double & max_y, const double &max_z);
	void setCellSize(const double cellSize);
	
	///////////////////////////////////////////////////////////////////
	//拿到叶子节点的Boundary，设置叶子节点的信息,得到叶子节点的数量
	void getLeafboundaryAndLeafNode(bool isTarget);


	///////////////////////////////////////////////////////////////////
	/**在doubleBuffer中导入指定的一帧*/
	void buildDBufferOctree(bool isTarget, ObjMesh * objMesh);
	void buildDBufferOctree(bool isTarget, vector<Vertex> &vertexList);
	
	/**根据给定两帧内容，创建双Buffer的八叉树*/
	void buildDBufferOctree(ObjMesh * referenceMesh, ObjMesh * targetMesh);
	void buildDBufferOctree(vector<Vertex> & referenceMeshVerList, vector<Vertex> & targetMeshVerList);

	//如果当前以及有两帧，再读入下一帧，然后更新八叉树
	void getNextTwoFrame(ObjMesh * newObjMesh);
	void getNextTwoFrame(vector<Vertex> &newVertexList);

	/**删除doubel buffer octree中指定的那一帧数据*/
	void clearOneFrame(bool isTarget);

	/**利用宽度优先遍历，得到byteSteam*/
	bool getByteSteamOfOctreeXOR(vector<char> & byteList_out);
	bool getByteSteamOfOctreeSigleFrame(vector<char> & byteList_out, bool isTarget);

	void clearDBufferOctree();

	/////////////////////////////
	//根据ByteList 以及swap的八叉树结构，恢复出相应reference 的八叉树结构
	/////////////////////////////
	
	/**根据byteStream来直接重构单帧的八叉树结构,
	* 前提是空的八叉树,该函数先清空八叉树
	*/
	bool recoverySigleFrameInDBufferOctree(vector<char> & byteListForSigleFrame, bool isTarget = false);

	/**根据byteStream来直接在已经建好第一帧的情况下，根据XOR byteStream来重构第二帧
	* 前提是已经重构好一帧了，通常是refe frame已经在double buffer octree中了
	*/
	bool recoverySecondFrameInDBufferOctreeXOR(vector<char> & byteStreamXOR, bool isTarget = true);


	/**根据swapFrameVerList作为reference frame 的点的信息，byteList是byteStream的信息，恢复出target frame 的八叉树的结果
	* 前提：空的double buffer octree
	* 方法: 根据swapFrameVerList_in，重建出单帧的八叉树结构，然后根据ByteList恢复出另外一个八叉树
	* */
	bool recoveryDBufferOctreeForTarget(vector<Vertex> &swapFrameVerList_in, vector<char> & byteList_in);

	bool getTargetOctreeByClearReference();
	
	

};
