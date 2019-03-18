#pragma once

#include <vector>
#include <set>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <Eigen/Dense>

#include "zjw_macro.h"

#ifdef SGWT_DEBUG
#include "zjw_sgwt_utils.h"
#endif //SGWT_DEBUG

#include "zjw_octree.h"
#include "util\zjw_math.h"
#include "util\zjw_obj.h"
#include "zjw_timer.h"
#include "zjw_dsaupd.h"
#include "zjw_kmeans.h"
#include "signalType.h"

using namespace Eigen;
using namespace std;

#ifdef USE_SPARSE
#include <eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
#endif // USE_SPARSE

//------------------------------ Class that holds your data.-----------------------

//节点包含点的数据
class Node
{
public:
	//每个在空间有一堆的 particle的位置信息
	vector<Vec3> pointPosList;

	//每个在空间有一堆的 particle的位置信息
	vector<Vec3> colorList;

	//每个在空间有一堆的 particle在obj中的序号
	vector<int> pointIdxList;

	//保存每个叶子节点下面8个象限的叶子节点的位置信息的指针
	vector<vector<Vec3*>> leafNodePos8Areas;
	//保存每个叶子节点下面8个象限的叶子节点的颜色信息的指针
	vector<vector<Vec3*>> leafNodeColor8Areas;

	//保存该叶子节点的在每个象限中的信号(true),如果没有信号，false signal (0,0,0)
	vector <bool> pos8Flag;

	//8个信号，每个象限一个。 每个node上的顶点位置的的信号。（x,y,z）
	vector<Vec3> pos8AreasSignal;

	//8个信号，每个象限一个。 每个node上的顶点位置的的信号。（x,y,z）
	vector<Vec3> color8AreasSignal;
};

//-----------------判断这个点在哪个叶子节点中-----------------
class AddPoints : public Octree<Node>::Callback
{
public:
	//保存当前这个点所在node的序号
	int nodeIdx;
	//需要被判断的点
	Vec3 point;
public:
	AddPoints(Vec3 point);
	~AddPoints();

	virtual bool operator()(const Vec3 min, const Vec3 max, Node& nodeData) {
		return true;
	}
	//返回false,终止递归遍历;返回true,继续递归遍历子节点
	virtual bool operator()(const Vec3 min, const Vec3 max, Octree<Node>::OctreeNode* currNode);
};

//-----------------遍历整个八叉树，把叶子节点的边界值保存到list中，并给叶子节点编号-----------------
class CallTraverseGetInfoSetLeaf : public Octree<Node>::Callback
{
public:
	//保存叶子节点的boundary ，用于渲染所用
	vector<Vec3> minVList;
	vector<Vec3> maxVList;
	//叶子节点的中间位置，用于kmeans
	vector<Vec3> *midVList;

	//保存Node上颜色的均值
	vector<Vec3> nodeColorList;
	//保存叶子节点的list,并给叶子节点编号
	vector<Octree<Node>::OctreeNode*> nodeList;
	//叶子节点的编号器
	int leafIncr;

public:
	CallTraverseGetInfoSetLeaf();
	~CallTraverseGetInfoSetLeaf();

	virtual bool operator()(const Vec3 min, const Vec3 max, Node& nodeData) {
		return true;
	}
	//给叶子节点编号
	virtual bool operator()(const Vec3 min, const Vec3 max, Octree<Node>::OctreeNode* currNode);
};

//----------------------把点云得到转换为矩阵----------------------------------

class CallTGetGraph : public Octree<Node>::Callback
{
public:
	//得到该节点周围邻居，及其权重，并保存到weight矩阵中。
	int leafIdx;
	Vec3 leafMidPoint;
	//不应该存储的预先设置的值，而应该是实际cell的宽度。
	Vec3 octCellSize;
#ifdef USE_EIGEN
	MatrixXd  * dMatPtr;
	MatrixXd  * weightAMatPtr;
#endif // USE_EIGEN

#ifdef USE_SPARSE
	SpMat * spLap;
	std::vector<T> * coeff;
#endif // USE_SPARSE

#ifdef  USE_ARPACK
	Dsaupd * laplacianMat;
#endif //  use_arpack

public:
#ifdef USE_EIGEN
	CallTGetGraph(MatrixXd  * dMatPtr, MatrixXd  * weightAMatPtr)
	{
		this->dMatPtr = dMatPtr;
		this->weightAMatPtr = weightAMatPtr;
	}
#endif // USE_EIGEN

#ifdef USE_SPARSE
	CallTGetGraph(std::vector<T> * coefficent, SpMat * sp = nullptr)
	{
		this->coeff = coefficent;
		this->spLap = sp;
	}
#endif // USE_SPARSE

#ifdef  USE_ARPACK
	CallTGetGraph(Dsaupd *laplacianMat)
	{
		this->laplacianMat = laplacianMat;
	}
#endif //  use_arpack

	~CallTGetGraph();

	void initParam(Vec3 & octreeCellSize);

	virtual bool operator()(const Vec3 min, const Vec3 max, Node& nodeData) {
		return true;
	}
	//得到相邻节点，及其权重
	virtual bool operator()(const Vec3 min, const Vec3 max, Octree<Node>::OctreeNode* currNode);
};

//-----------------------------------------------------
class PcsOctree
{
public:
	Vec3 minPos;
	Vec3 maxPos;
	Vec3 cellSize;

	//ObjMesh * objMesh;
	Octree<Node> *pcsOct;
	CallTraverseGetInfoSetLeaf * ctLeaf;
	CallTGetGraph *ctGraph;
	//对八叉树的叶子节点进行k-means聚类
	KMeans * kmeans;

	//sgwt 中的参数
	//Chebyshev的系数，近似的M阶的相数
	int m;
	//sclae的个数，不包括h
	int Nscales;

	//邻接权重矩阵，以及D矩阵，求解特征值和特征向量的存储矩阵。
	int nodeNum;

	//8个象限中的信号
	vector<VectorXd> posSignalX;
	vector<VectorXd> posSignalY;
	vector<VectorXd> posSignalZ;

#ifdef USE_EIGEN
	MatrixXd  dMat;
	MatrixXd  weightAMat;
	MatrixXd  LaplacianMat;
	MatrixXd  eigenVecMat;
	MatrixXd  eigenValMat;
#endif // use_ei

#ifdef USE_SPARSE
	SpMat * spLaplacian;
	std::vector<T> * coefficients;
	//sgwt 的对象
#ifdef SGWT_DEBUG
	SgwtCheby *fastSgwt;
#endif //SGWT_DEBUG
#endif // USE_SPARSE

#ifdef USE_ARPACK
	Dsaupd laplacianMat;
#endif // use_arpa

public:
	PcsOctree(int _m, int _Nscales);
	~PcsOctree();
	//设置参数，重新初始化八叉树
	void setParam(Vec3 min, Vec3 max, Vec3 cellSize);
	//建立八叉树
	void buildPcsOctFrmPC(ObjMesh * objeMesh);
	//拿到叶子节点的Boundary，设置叶子节点的信息
	void getLeafboundary();

#ifdef RELATIVE_DIS_SIGNAL
	void changeNodePointsToRelativeDis();
#endif // RELATIVE_DIS_SIGNAL

	//把八叉树node中的关于这个node所有点的坐标改成相对坐标。

	//拿到graph的两个矩阵
	//初始化矩阵等信息，在getGraphMat之前调用
	void initMat();
	//得到矩阵
	void getGraphMat();

	//得到矩阵的特征向量和特征值
	void getMatEigenVerValue();

	//清空八叉树
	void clearOct();

	//判断顶点point 在包围盒中的8个象限中的哪一个
	int judege8Aeros(Vec3 &min, Vec3 &max, Vec3 &point);
	int judege8Aeros(Vec3& mid, Vec3 &point);

	//把每个叶子节点中的点划分到8个子象限中，并计算每个叶子节点的信号，保存在nodeData中
	void setPointTo8Areas();
	//计算叶子节点在每个象限上的信号
	void getLeafSignal();

	//拿到叶子节点上的平均颜色
	vector<Vec3> & getLeafAvgColor() const;

	//把8个象限的信号（x,y,z,r,g,b）放到向量中:并封装到vector中
	//8个vectorXd， 分别表示不同象限的信号。每个vector表示的是：所有在节点在这个象限中的该型号的具体的值
	vector<VectorXd> getSignalF(SignalType sType);

	//调用getSignalF得到所有信号，利用sgwt求出所有的系数，并保存在sgwt中。
	bool getAllSignalAndSaveSGWTCoeff();

	//返回这个节点在特定信号下面所有象限信号，返回值应该是8个系数的向量，对应每个象限
	bool getSignalF(SignalType sType, int nodeIdx, VectorXd * sigVQ_out);

	//废弃,输入信号的值和象限的类型(0-7),返回node上在所有持尺度上的系数
	vector<VectorXd> getSgwtCoeffWS(SignalType type, int quadrant);

	//得到这个node的feature vector 维数 8*6*5
	bool getFeatureVector(int nodeIdx, VectorXd* featureVector);
	//加速
	bool getFeatureVector2(int nodeIdx, VectorXd* featureVector);

	//废弃，x,y，z所有的信号
	void getSgwtCoeffWS();

	//kmeans
	void doKmeans(int clusterNum);
	//判断当前点在哪个叶子结点上,返回的是叶子节点的序号
	int judegePointToLeafNode(Vec3* point, int& idx);


	//拿到nodeidx这个点的two hop。拿到的点不包括自己，因为是无向图，所以这里的two hop包括一步的
	void getTwoHopNeighborhood(int nodeIdx, set<int> * nodeList_out, SpMat * spLaplacian);

	//two hop包括一步的，上面的接口是不包括一步的。
	void getTwoHopNeighborhoodWithOneStep(int nodeIdx, set<int> * nodeList_out, SpMat * spLaplacian);
	
	//对外的接口
	//void
	//=============== test ====================
	void printMat();
};
