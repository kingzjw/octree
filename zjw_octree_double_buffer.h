#pragma once
/*!
* \file zjw_octree.h
*
* \author King
* \date 7月 2018
*  this head file and the implements c++ souce code has set the basis octree build and opration
*  interface. which implements in real application should look at zjw_pcs_octree.h
*/

#include <assert.h>
#include <vector>
#include "util/zjw_macro.h"
#include <util/zjw_math.h>

/*
* 作用：
*	定义宏，得到点所在象限位置
* 
* 输入：
*	bit: 表示对x,y,z中的哪一维度进行判断
*   p: 表示的当前需要进行判断的点云帧中的顶点
*   mid: 表示的是当前Node的中点的坐标，需要根据这个坐标分出8个象限
* 输出：
*	i: 表示当前p顶点，应该在当前的node上的哪一个象限
*   newMin: 表示新的象限边界的最小点的值
*   newMax: 表示新的象限边界的最大点的值
*/

#define COMPUTE_QUADRANT(i, bit, p, mid, newMin, newMax) \
if (p >= mid)         \
{                     \
    i |= bit;         \
    newMin = mid;     \
}                     \
else                  \
{                     \
    newMax = mid;     \
}

////////////////////////////////////////////////////////////////////////////////////////////////
//double buffer的八叉树的结点
////////////////////////////////////////////////////////////////////////////////////////////////

template <class NodeDataType>
class DoubleBufferOctree;

template <class NodeDataType>
class OctreeDoubelBufferNode
{
private:
	
	//存放reference target frame的指针
	DoubleBufferOctree<NodeDataType> * oct;

	//非叶子节点是-1，其他的从0,1,2,3一次编号.这个参数要被使用，需要现给它值。
	int leafFlag;
	
	//当前节点的边界值的最小值和最大值
	Vec3 min;
	Vec3 max;


public:
	//元素是指针
	OctreeDoubelBufferNode* children[8];
	//这个Node上的相关的点云数据
	NodeDataType * nodeData[2];
	//标记
	bool flag[2][8];

public:
	OctreeDoubelBufferNode(DoubleBufferOctree<NodeDataType> * oct);

	/** isTarget: 是什么node类型，然后需要对相应的nodeData 创建对象	*/
	OctreeDoubelBufferNode(DoubleBufferOctree<NodeDataType> * oct, bool isTarget);
	virtual ~OctreeDoubelBufferNode();

	NodeDataType* getTargetFrameNodeDataPtr();
	NodeDataType* getRefrenceFrameNodeDataPtr();

	//设置当前节点，有子节点（true）
	int getTargetFlagIdx();
	int getReferenceFlagIdx();
	void setFlag(bool isTarget, int index);
	//设置当前节点，没有子节点（false）
	void resetFlag(bool isTarget, int index);
	//初始化指定的Flag
	void initFlag(bool isTarget);
	//初始化所有的falg
	void initFlag();

	//初始化 特定帧的 nodedata
	NodeDataType * getNodeData(bool isTarget);
	void setNodeData(bool isTarget);
	void resetNodeData(bool isTarget);
	//检测nodeData是否创建对象 true:表示Null
	bool checkNodeDataNull(bool isTarget);

	//得到孩子的指针,如果是空指针，一定要特别注意，不能把空指针赋值给一个空指针//
	OctreeDoubelBufferNode* getChildren(int index);
	Vec3 getMinPos();
	Vec3 getMaxPos();
	void setMinPos(Vec3 minP);
	void setMaxPos(Vec3 maxP);
	void setLeafFlag(int i);

	//判断这个node是否是叶子节点
	bool judegeLeafNode();

	//XOR根据flag中的内容，计算出这个byte
	char getChildByteXOR();
	char getChildByteSigleFrame(bool isTarget);

	//清理这一帧的相关数据
	void resetFrame(bool isTarget);
};

////////////////////////////////////////////////////////////////////////////////////////////////
//double buffer 的 八叉树
////////////////////////////////////////////////////////////////////////////////////////////////

template <class NodeDataType>
class DoubleBufferOctree
{
public:
	friend class OctreeDoubelBufferNode<NodeDataType>;

	//min的顶点和Max的顶点确定 节点的空间范围
	Vec3 min;
	Vec3 max;
	// cellSize is the minimum subdivided cell size.
	Vec3 cellSize;
	//根节点
	OctreeDoubelBufferNode<NodeDataType>* root;

	//双Buffer的序号标志
	int targetFrameId;
	int refrenceFrameId;

public:

	DoubleBufferOctree(Vec3 min, Vec3 max, Vec3 cellSize,int targetFrameId = 1 , int refrenceFrameId = 0) ;
	virtual ~DoubleBufferOctree();

	void setTargetFrameId(int id);
	void setReferenceFrameId(int id);
	
	//创建指定父亲节点的还是节点，并返回这个孩子节点的指针
	OctreeDoubelBufferNode<NodeDataType>* createChildNode(bool isTarget,OctreeDoubelBufferNode<NodeDataType>* parentNode, int childIdx);
	
	//decode: 根据字节，恢复这个node的孩子信息。前提是db octree是空的，没有任何信息
	bool recoveryDBufferOctreeNodeForSigleFrame(OctreeDoubelBufferNode<NodeDataType>* parentNode, const char byte, bool isTarget =false);

	//decode: 根据字节，恢复这个node的孩子信息。前提是已经有单帧存在这个八叉树上了，用于恢复另外一帧
	bool recoveryDBufferOctreeNodeXOR(OctreeDoubelBufferNode<NodeDataType>* parentNode, const char byteXOR);

	
	//////////////////////////////////////////////
	//遍历所有的Octree中做操作的一个类。
	//////////////////////////////////////////////

	class Callback
	{
	public:
		// Return value: true = continue; false = abort.
		virtual bool operator()(const Vec3 min, const Vec3 max, NodeDataType& nodeData) { return true; }

		// Return value: true = continue; false = abort.
		virtual bool operator()(const Vec3 min, const Vec3 max, OctreeDoubelBufferNode<NodeDataType>* currNode) { return true; }
	};

	//根据深度遍历的方法，把特定点pos插入到八叉树的叶子节点中，直到叶子节点的size小于 cellSize.
	/** 因为是double buffer，所以需要对target 或者是reference 进行指定
	* 返回的是叶子节点中的Node，可以进行叶子节点的处理
	* */
	NodeDataType* getCell(bool isTarget, const Vec3 ppos, Callback* callback = NULL);

	//深度优先递归遍历整个树.并进行某种操作。（先操作后深入，深入是有条件判断的）
	void traverse(Callback* callback);
	//深度优先遍历整个树，所有的节点都会遍历到，从下往上进行操作。	
	void traverseDepthFirst(Callback* callback);
	//广度优先遍历所有节点,得到所有列表的参数
	void traverseBreathFirst(std::vector< OctreeDoubelBufferNode<NodeDataType>*> &allNodeList_out);
	void traverseBreathFirstSigleFrame(std::vector< OctreeDoubelBufferNode<NodeDataType>*> &allNodeList_out,bool isTarget_in);

	//删除整个树，释放所有节点的空间
	void clear();
	

protected:

	//和traverse函数配套使用
	void depthFirstTraverseRecursive(Callback* callback, const Vec3& currMin, const Vec3& currMax, OctreeDoubelBufferNode<NodeDataType>* currNode);
	
	//和traverseDepthFirst配套使用
	void traverseDepthFirstRecursive(Callback* callback, const Vec3& currMin, const Vec3& currMax, OctreeDoubelBufferNode<NodeDataType>* currNode);

};


/////////////////////////////////////////////////////////////////////////////////////////


template<class NodeDataType>
OctreeDoubelBufferNode<NodeDataType>::OctreeDoubelBufferNode(DoubleBufferOctree<NodeDataType> * oct)
{
	this->oct = oct;

	leafFlag = -1;
	//初始化两个node data
	nodeData[0] = 0;
	nodeData[1] = 0;

	//初始化八个孩子的结点
	for (int i = 0; i < 8; i++)
		children[i] = 0;

	initFlag();

}

template<class NodeDataType>
OctreeDoubelBufferNode<NodeDataType>::OctreeDoubelBufferNode(DoubleBufferOctree<NodeDataType> * oct, bool isTarget)
{
	this->oct = oct;

	leafFlag = -1;
	
	//初始化两个node data
	nodeData[0] = 0;
	nodeData[1] = 0;

	//初始化八个孩子的结点
	for (int i = 0; i < 8; i++)
		children[i] = 0;

	initFlag();
	setNodeData(isTarget);
}

template<class NodeDataType>
inline OctreeDoubelBufferNode<NodeDataType>::~OctreeDoubelBufferNode()
{
	for (int i = 0; i < 8; i++)
	{
		if (children[i])
			delete children[i];
	}

	if (nodeData[0])
		delete nodeData[0];

	if (nodeData[1])
		delete nodeData[1];
}

template<class NodeDataType>
void OctreeDoubelBufferNode<NodeDataType>::initFlag()
{
	//init flag
	for (int i = 0; i < 8; i++)
	{
		flag[0][i] = 0;
		flag[1][i] = 0;
	}
}

template<class NodeDataType>
void OctreeDoubelBufferNode<NodeDataType>::initFlag(bool isTarget)
{
	if (isTarget)
	{
		for (int i = 0; i < 8; i++)
		{
			flag[oct->targetFrameId][i] = 0;
		}
	}
	else
	{
		for (int i = 0; i < 8; i++)
		{
			flag[oct->refrenceFrameId][i] = 0;
		}
	}
}

template<class NodeDataType>
NodeDataType * OctreeDoubelBufferNode<NodeDataType>::getNodeData(bool isTarget)
{
	if (isTarget)
		return nodeData[oct->targetFrameId];
	else
		return nodeData[oct->refrenceFrameId];
}

template<class NodeDataType>
void OctreeDoubelBufferNode<NodeDataType>::setNodeData(bool isTarget)
{
	if (isTarget)
		nodeData[oct->targetFrameId] = new NodeDataType();
	else
		nodeData[oct->refrenceFrameId] = new NodeDataType();
}

template<class NodeDataType>
void OctreeDoubelBufferNode<NodeDataType>::resetNodeData(bool isTarget)
{
	if (isTarget)
	{
		if (nodeData[oct->targetFrameId])
			delete nodeData[oct->targetFrameId];
	}
	else
	{
		if (nodeData[oct->refrenceFrameId])
			delete nodeData[oct->refrenceFrameId];
	}
}

template<class NodeDataType>
bool OctreeDoubelBufferNode<NodeDataType>::checkNodeDataNull(bool isTarget)
{
	if (isTarget)
	{
		if (nodeData[oct->targetFrameId])
			return false;
	}
	else
	{
		//test
		//cout << oct << endl;
		//cout << oct->refrenceFrameId << endl;
		//end test
		if (nodeData[oct->refrenceFrameId])
			return false;
	}
	return true;
}

template<class NodeDataType>
OctreeDoubelBufferNode<NodeDataType> * OctreeDoubelBufferNode<NodeDataType>::getChildren(int index)
{
	return children[index];
}

template<class NodeDataType>
Vec3 OctreeDoubelBufferNode<NodeDataType>::getMinPos()
{
	return this->min;
}

template<class NodeDataType>
Vec3 OctreeDoubelBufferNode<NodeDataType>::getMaxPos()
{
	return this->max;
}

template<class NodeDataType>
inline void OctreeDoubelBufferNode<NodeDataType>::setMinPos(Vec3 minP)
{
	this->min = minP;
}

template<class NodeDataType>
inline void OctreeDoubelBufferNode<NodeDataType>::setMaxPos(Vec3 maxP)
{
	this->max = maxP;
}

template<class NodeDataType>
void OctreeDoubelBufferNode<NodeDataType>::setLeafFlag(int i)
{
	leafFlag = i;
}

template<class NodeDataType>
bool OctreeDoubelBufferNode<NodeDataType>::judegeLeafNode()
{
	if (this->getChildren(0) || this->getChildren(1) || this->getChildren(2) || this->getChildren(3)
		|| this->getChildren(4) || this->getChildren(5) || this->getChildren(6) || this->getChildren(7))
	{
		return false;
	}
	return true;
}

template<class NodeDataType>
char OctreeDoubelBufferNode<NodeDataType>::getChildByteXOR()
{
	bitset<8> b;

	for (int i = 0; i < 8; i++)
	{
		if (flag[0][i] == flag[1][i])
		{
			b.set(i, 0);
		}
		else
		{
			b.set(i, 1);
		}
	}
	char ch = int(b.to_ulong());
	return ch;
}

template<class NodeDataType>
inline char OctreeDoubelBufferNode<NodeDataType>::getChildByteSigleFrame(bool isTarget)
{
	int  idx = -1;
	if (isTarget)
	{
		idx = getTargetFlagIdx();
	}
	else
	{
		idx = getReferenceFlagIdx();
	}

	bitset<8> b;
	b.reset();
	for (int i = 0; i < 8; i++)
	{
		if (flag[idx][i])
		{
			b.set(i, 1);
		}
	}
	return int(b.to_ulong());
}

template<class NodeDataType>
inline NodeDataType * OctreeDoubelBufferNode<NodeDataType>::getTargetFrameNodeDataPtr()
{
	if (oct->targetFrameId == 0)
		return nodeData[0];
	else if (oct->targetFrameId == 1)
		return nodeData[1];
}

template<class NodeDataType>
inline NodeDataType * OctreeDoubelBufferNode<NodeDataType>::getRefrenceFrameNodeDataPtr()
{
	if (oct->refrenceFrameId == 0)
		return nodeData[0];
	else if (oct->refrenceFrameId == 1)
		return nodeData[1];
}

template<class NodeDataType>
int OctreeDoubelBufferNode<NodeDataType>::getTargetFlagIdx()
{
	return oct->targetFrameId;
}

template<class NodeDataType>
int OctreeDoubelBufferNode<NodeDataType>::getReferenceFlagIdx()
{
	return oct->refrenceFrameId;
}

template<class NodeDataType>
void OctreeDoubelBufferNode<NodeDataType>::setFlag(bool isTarget, int index)
{
	if (isTarget)
		flag[oct->targetFrameId][index] = true;
	else
		flag[oct->refrenceFrameId][index] = true;
}

template<class NodeDataType>
void OctreeDoubelBufferNode<NodeDataType>::resetFlag(bool isTarget, int index)
{
	if (isTarget)
		flag[oct->targetFrameId][index] = false;
	else
		flag[oct->refrenceFrameId][index] = false;
}

template<class NodeDataType>
void OctreeDoubelBufferNode<NodeDataType>::resetFrame(bool isTarget)
{
	initFlag(isTarget);
	resetNodeData(isTarget);
}

/////////////////////////////////////////////////////////////////////////////////////////

template<class NodeDataType>
inline DoubleBufferOctree<NodeDataType>::DoubleBufferOctree(Vec3 min, Vec3 max, Vec3 cellSize, int targetFrameId, int refrenceFrameId)
{
	assert(targetFrameId == 0 || targetFrameId == 1);
	assert(refrenceFrameId == 0 || refrenceFrameId == 1);

	this->min = min;
	this->max = max;
	this->cellSize = cellSize;
	this->targetFrameId = targetFrameId;
	this->refrenceFrameId = refrenceFrameId;

	root = new OctreeDoubelBufferNode<NodeDataType>(this);
	//min的顶点和Max的顶点确定 节点的空间范围
	root->setMinPos(min);
	root->setMaxPos(max);
}

template<class NodeDataType>
inline DoubleBufferOctree<NodeDataType>::~DoubleBufferOctree()
{
	if (root)
		delete root;
}

template<class NodeDataType>
void DoubleBufferOctree<NodeDataType>::setTargetFrameId(int id)
{
	assert(id == 0 || id == 1);
	//一个是0，一个是1
	targetFrameId = id;
	refrenceFrameId = (id + 1) % 2;
}

template<class NodeDataType>
void DoubleBufferOctree<NodeDataType>::setReferenceFrameId(int id)
{
	assert(id == 0 || id == 1);
	//一个是0，一个是1
	targetFrameId = (id + 1) % 2;
	refrenceFrameId = id;
}

template<class NodeDataType>
OctreeDoubelBufferNode<NodeDataType>* DoubleBufferOctree<NodeDataType>::createChildNode(bool isTarget,
	OctreeDoubelBufferNode<NodeDataType>* parentNode, int childIdx)
{
	assert(parentNode);
	assert(childIdx > -1 && childIdx < 8);
	
	if (!parentNode->getChildren(childIdx))
	{
		//没有这个孩子节点，那么创建这个节点
		parentNode->children[childIdx] = new OctreeDoubelBufferNode<NodeDataType>(this, isTarget);
		OctreeDoubelBufferNode<NodeDataType> * temp = parentNode->getChildren(childIdx);

		//根据父亲节点的边界，得到孩子节点的边界
		Vec3 newMin = parentNode->getMinPos();
		Vec3 newMax = parentNode->getMaxPos();
		Vec3 mid = (parentNode->getMinPos() + parentNode->getMaxPos())/2;
		switch (childIdx)
		{
		case 0:
			newMax = mid;
			break;
		case 1:
			newMin.x = mid.x;
			newMax.y = mid.y;
			newMax.z = mid.z;
			break;
		case 2:
			newMin.y = mid.y;
			newMax.x = mid.x;
			newMax.z = mid.z;
			break;
		case 3:
			newMin.x = mid.x;
			newMin.y = mid.y;
			newMax.z = mid.z;
			break;
		case 4:
			newMin.z = mid.z;
			newMax.x = mid.x;
			newMax.y = mid.y;
			break;
		case 5:
			newMin.x = mid.x;
			newMin.z = mid.z;
			newMax.y = mid.y;
			break;
		case 6:
			newMin.y = mid.y;
			newMin.z = mid.z;
			newMax.x = mid.x;
			break;
		case 7:
			newMin = mid;
			break;
		default:
			cout << "do not support the child idx!" << endl;
			break;
		}

		temp->setMinPos(newMin);
		temp->setMaxPos(newMax);
		//flag
		if (isTarget)
			parentNode->flag[targetFrameId][childIdx] = true;
		else
			parentNode->flag[refrenceFrameId][childIdx] = true;

		//创建该叶子节点的nodedata
		temp->setNodeData(isTarget);


		return temp;
	}
	else
	{
		//已经有这个孩子节点，那么创建的是另一帧的情况

		//flag
		if (isTarget)
		{
			parentNode->flag[targetFrameId][childIdx] = true;
		}
		else
		{
			parentNode->flag[refrenceFrameId][childIdx] = true;
		}
		//node data： 叶子节点存在，需要创建这个叶子节点上，针对该帧的nodedata
		parentNode->getChildren(childIdx)->setNodeData(isTarget);

		//if (parentNode->checkNodeDataNull(isTarget))
			//parentNode->setNodeData(isTarget);

		return parentNode->getChildren(childIdx);
	}
}

template<class NodeDataType>
inline bool DoubleBufferOctree<NodeDataType>::recoveryDBufferOctreeNodeForSigleFrame(OctreeDoubelBufferNode<NodeDataType>* parentNode, const char byte, bool isTarget)
{
	assert(parentNode);

	bitset<8> bits = byte;
	//处理当前节点的的8个孩子
	for (int b_it=0; b_it < 8; b_it++)
	{
		//这一位是1，创建这个节点
		if (bits[b_it])
		{
			//创建这个节点的，第b-it个孩子节点
			createChildNode(isTarget, parentNode, b_it);
		}
	}
	return true;
}

template<class NodeDataType>
bool DoubleBufferOctree<NodeDataType>::recoveryDBufferOctreeNodeXOR(OctreeDoubelBufferNode<NodeDataType>* parentNode, const char byte)
{
	assert(parentNode);

	bitset<8> bits = byte;
	//处理当前节点的的8个孩子
	for (int b_it = 0; b_it < 8; b_it++)
	{
		//这一位是1，那么说明孩子节点上是不一样，需要进行处理
		if (bits[b_it])
		{
			if (!parentNode->flag[refrenceFrameId][b_it])
			{
				//swapframe 在这个节点上是没有孩子的，targetFrame是有孩子不用处理
				createChildNode(true, parentNode, b_it);
			}
			else
			{
				//swapframe 在这个节点上是有孩子的，targetFrame是没有孩子，不用处理
			}
		}
		else
		{
			//这一位是0，说明孩子节点是一样的

			if (parentNode->flag[refrenceFrameId][b_it])
			{
				//swapframe 在这个节点上是有孩子的，targetFrame是有孩子，用处理
				createChildNode(true, parentNode, b_it);
			}
			else
			{
				//swapframe 在这个节点上是没有孩子的，不用处理

			}
		}
	}
	return true;
}

template<class NodeDataType>
NodeDataType * DoubleBufferOctree<NodeDataType>::getCell(bool isTarget, const Vec3 ppos, Callback * callback)
{
	assert(ppos >= min && ppos <= max);

	Vec3 currMin(min);
	Vec3 currMax(max);
	Vec3 delta = max - min;

	if (!root)
	{
		//root是一早就被初始化了
		root = new OctreeDoubelBufferNode<NodeDataType>(this, isTarget);
	}
	//如果root已经分配，那么检查特定frame的分配情况.因为有两个frame需要进行管理
	if (root->checkNodeDataNull(isTarget))
		root->setNodeData(isTarget);

	OctreeDoubelBufferNode<NodeDataType> * currNode = root;

	//判断是否需要继续往八叉树里面走
	while (delta >= cellSize)
	{
		bool shouldContinue = true;

		if (callback)
			shouldContinue = callback->operator()(currMin, currMax, *(currNode->getNodeData(isTarget)));
		if (!shouldContinue)
			break;

		//得到当前空间的中心点
		Vec3 mid = (delta * 0.5f) + currMin;

		//初始化下一个子空间的最小位置和最大位置。
		Vec3 newMin(currMin);
		Vec3 newMax(currMax);
		//index记录点所在象限位置。
		int index = 0;

		//判断ppos这个点，在当前的这个node范围中，应该属于哪个象限,并得到子空间的位置
		COMPUTE_QUADRANT(index, 1, ppos.x, mid.x, newMin.x, newMax.x);
		COMPUTE_QUADRANT(index, 2, ppos.y, mid.y, newMin.y, newMax.y);
		COMPUTE_QUADRANT(index, 4, ppos.z, mid.z, newMin.z, newMax.z);

		//如果这个子节点都没有分配，那么进行分配
		
		if (!(currNode->getChildren(index)))
		{
			currNode->children[index] = new OctreeDoubelBufferNode<NodeDataType>(this, isTarget);
			//设置flag
			//currNode->setFlag(isTarget, index);
		}
		//如果子空间节点已经分配，那么检查特定frame中的子节点的nodedata的分配情况
		if (currNode->children[index]->checkNodeDataNull(isTarget))
			currNode->children[index]->setNodeData(isTarget);

		currNode->setFlag(isTarget, index);

		//currNode = currNode->getChildren(index);
		currNode = currNode->children[index];
		//每个节点保留自己的范围
		currNode->setMinPos(newMin);
		currNode->setMaxPos(newMax);

		currMin = newMin;
		currMax = newMax;
		delta = currMax - currMin;
	}
	return currNode->getNodeData(isTarget);
}

template<class NodeDataType>
void DoubleBufferOctree<NodeDataType>::traverse(Callback * callback)
{
	assert(callback);
	depthFirstTraverseRecursive(callback, min, max, root);
}

template<class NodeDataType>
void DoubleBufferOctree<NodeDataType>::traverseDepthFirst(Callback * callback)
{
	assert(callback);
	traverseDepthFirstRecursive(callback, min, max, root);
}

//allNodeList_out 返回的值中包含叶子节点。
template<class NodeDataType>
void DoubleBufferOctree<NodeDataType>::traverseBreathFirst(std::vector<OctreeDoubelBufferNode<NodeDataType>*> &allNodeList_out)
{

	//宽度优先遍历，得到所有节点的一个序列
	allNodeList_out.push_back(root);
	for (int i = 0; i < allNodeList_out.size(); i++)
	{
		for (int child_it = 0; child_it < 8; child_it++)
		{
			if (allNodeList_out[i]->getChildren(child_it))
			{
				allNodeList_out.push_back(allNodeList_out[i]->getChildren(child_it));
			}
		}
	}
}

template<class NodeDataType>
inline void DoubleBufferOctree<NodeDataType>::traverseBreathFirstSigleFrame(std::vector<OctreeDoubelBufferNode<NodeDataType>*>& allNodeList_out, bool isTarget_in)
{
	//宽度优先遍历，得到所有节点的一个序列
	if (isTarget_in)
	{
		allNodeList_out.push_back(root);
		for (int i = 0; i < allNodeList_out.size(); i++)
		{
			for (int child_it = 0; child_it < 8; child_it++)
			{
				//target 上的叶子节点是OK的
				if (allNodeList_out[i]->flag[allNodeList_out[i]->getTargetFlagIdx()][child_it])
				{
					allNodeList_out.push_back(allNodeList_out[i]->getChildren(child_it));
				}
			}
		}
	}
	else
	{
		allNodeList_out.push_back(root);
		for (int i = 0; i < allNodeList_out.size(); i++)
		{
			for (int child_it = 0; child_it < 8; child_it++)
			{
				if (allNodeList_out[i]->flag[allNodeList_out[i]->getReferenceFlagIdx()][child_it])
				{
					allNodeList_out.push_back(allNodeList_out[i]->getChildren(child_it));
				}
			}
		}
	}
}

template<class NodeDataType>
void DoubleBufferOctree<NodeDataType>::clear()
{
	//为空返回
	if (!currNode)
		return;
	//判断是否符合其他终止条件
	bool shouldContinue = callback->operator()(currMin, currMax, currNode);
	if (!shouldContinue)
		return;

	//mid 得到cell中间的点
	Vec3 delta = currMax - currMin;
	Vec3 mid = (delta * 0.5f) + currMin;
	//分成子八个节点
	depthFirstTraverseRecursive(callback, currMin, mid, currNode->children[0]);
	depthFirstTraverseRecursive(callback, Vec3(mid.x, currMin.y, currMin.z), Vec3(currMax.x, mid.y, mid.z), currNode->children[1]);
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, mid.y, currMin.z), Vec3(mid.x, currMax.y, mid.z), currNode->children[2]);
	depthFirstTraverseRecursive(callback, Vec3(mid.x, mid.y, currMin.z), Vec3(currMax.x, currMax.y, mid.z), currNode->children[3]);
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, currMin.y, mid.z), Vec3(mid.x, mid.y, currMax.z), currNode->children[4]);
	depthFirstTraverseRecursive(callback, Vec3(mid.x, currMin.y, mid.z), Vec3(currMax.x, mid.y, currMax.z), currNode->children[5]);
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, mid.y, mid.z), Vec3(mid.x, currMax.y, currMax.z), currNode->children[6]);
	depthFirstTraverseRecursive(callback, mid, currMax, currNode->children[7]);
}

template<class NodeDataType>
void DoubleBufferOctree<NodeDataType>::depthFirstTraverseRecursive(Callback * callback, const Vec3 & currMin, const Vec3 & currMax, OctreeDoubelBufferNode<NodeDataType> * currNode)
{
	//为空返回
	if (!currNode)
		return;
	//判断是否符合其他终止条件
	bool shouldContinue = callback->operator()(currMin, currMax, currNode);
	if (!shouldContinue)
		return;

	//mid 得到cell中间的点
	Vec3 delta = currMax - currMin;
	Vec3 mid = (delta * 0.5f) + currMin;
	//分成子八个节点
	depthFirstTraverseRecursive(callback, currMin, mid, currNode->getChildren(0));
	depthFirstTraverseRecursive(callback, Vec3(mid.x, currMin.y, currMin.z), Vec3(currMax.x, mid.y, mid.z), currNode->getChildren(1));
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, mid.y, currMin.z), Vec3(mid.x, currMax.y, mid.z), currNode->getChildren(2));
	depthFirstTraverseRecursive(callback, Vec3(mid.x, mid.y, currMin.z), Vec3(currMax.x, currMax.y, mid.z), currNode->getChildren(3));
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, currMin.y, mid.z), Vec3(mid.x, mid.y, currMax.z), currNode->getChildren(4));
	depthFirstTraverseRecursive(callback, Vec3(mid.x, currMin.y, mid.z), Vec3(currMax.x, mid.y, currMax.z), currNode->getChildren(5));
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, mid.y, mid.z), Vec3(mid.x, currMax.y, currMax.z), currNode->getChildren(6));
	depthFirstTraverseRecursive(callback, mid, currMax, currNode->getChildren(7));
}

template<class NodeDataType>
void DoubleBufferOctree<NodeDataType>::traverseDepthFirstRecursive(Callback * callback, const Vec3 & currMin, const Vec3 & currMax, OctreeDoubelBufferNode<NodeDataType>* currNode)
{
	//为空返回
	if (!currNode)
		return;

	//mid 得到cell中间的点
	Vec3 delta = currMax - currMin;
	Vec3 mid = (delta * 0.5f) + currMin;
	//分成子八个节点
	depthFirstTraverseRecursive(callback, currMin, mid, currNode->getChildren(0));
	depthFirstTraverseRecursive(callback, Vec3(mid.x, currMin.y, currMin.z), Vec3(currMax.x, mid.y, mid.z), currNode->getChildren(1));
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, mid.y, currMin.z), Vec3(mid.x, currMax.y, mid.z), currNode->getChildren(2));
	depthFirstTraverseRecursive(callback, Vec3(mid.x, mid.y, currMin.z), Vec3(currMax.x, currMax.y, mid.z), currNode->getChildren(3));
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, currMin.y, mid.z), Vec3(mid.x, mid.y, currMax.z), currNode->getChildren(4));
	depthFirstTraverseRecursive(callback, Vec3(mid.x, currMin.y, mid.z), Vec3(currMax.x, mid.y, currMax.z), currNode->getChildren(5));
	depthFirstTraverseRecursive(callback, Vec3(currMin.x, mid.y, mid.z), Vec3(mid.x, currMax.y, currMax.z), currNode->getChildren(6));
	depthFirstTraverseRecursive(callback, mid, currMax, currNode->getChildren(7));

	//执行相应的操作代码块
	callback->operator()(currMin, currMax, currNode);

}


