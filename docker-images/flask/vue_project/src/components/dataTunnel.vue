<script setup>
import { computed } from 'vue';
import { useEdge, BezierEdge, MarkerType } from '@vue-flow/core';

const { edge } = useEdge();

edge.style = {
    stroke: 'black',
    'stroke-width': '4px'
};

edge.markerEnd = {
  type: MarkerType.Arrow,
  color: "black",
  strokeWidth: "2px" 
}

edge.animated = computed(() => {
    return edge.sourceNode.active;
});

</script>

<template>
    <BezierEdge/>
</template>