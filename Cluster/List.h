//����ڵ�
struct ListNode
{
	int m_nValue;
	ListNode* m_pNext;
};

//������β����ӽڵ�
void AddToTail(ListNode** pHead, int value);
//���������ҵ���һ������ĳֵ�Ľڵ㲢ɾ���Ľڵ�
void RemoveNode(ListNode** pHead, int value);

void PrintNode(ListNode *pHead);