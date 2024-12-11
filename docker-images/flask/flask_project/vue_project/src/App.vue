<script setup>
import { ref, markRaw } from 'vue';
import { VueFlow, useVueFlow } from '@vue-flow/core';
import { MiniMap } from '@vue-flow/minimap'
import '@vue-flow/minimap/dist/style.css'
import { Controls } from '@vue-flow/controls'
import '@vue-flow/controls/dist/style.css'
import { Background } from '@vue-flow/background'
import dataTank from './components/dataTank.vue';
import runJob from './components/runJob.vue';
import dataTunnel from './components/dataTunnel.vue';
import { validTargets, initRunJobNode, initDataTankNode } from './components/utils.js';

const nodeTypes = {
  dataTankNode: markRaw(dataTank),
  runJobNode: markRaw(runJob)
};

const edgeTypes = {
  dataTunnelEdge: markRaw(dataTunnel)
};

const nodes = ref([
  initDataTankNode(
    'target file', 650, 50,
    [ { type: 'file', name: 'target file' } ],
    'target reads aligned to references (.fq[.gz])'
  ),
  initDataTankNode(
    'pair file', 650, 150,
    [ { type: 'file', name: 'pair file' } ],
    'auxiliary reads for demultiplexing target reads (.fq[.gz])'
  ),
  initRunJobNode(
    'removeDuplicates', 1000, 125,
    'remove duplicated reads in target\\pair file'
  ),
  initDataTankNode(
    'file without duplicates', 1300, 100,
    [ { type: 'file', name: 'file without duplicates' } ],
    'file of side-by-side target and pair reads without duplicates'
  ),

  initDataTankNode(
    'target spliter', 650, 250,
    [ { type: 'file', name: 'target spliter' } ],
    'reads for demultiplexing target reads (.fa)'
  ),
  initRunJobNode(
    'buildSpliter/target spliter', 1000, 250,
    'build bowtie2 index for target spliter'
  ),
  initDataTankNode(
    'target spliter index', 1300, 200,
    [
      { type: 'value', name: 'minimal alignment score of target spliter', value: null },
      { type: 'files', name: 'target spliter index', value: 6 }
    ],
    'bowtie2 index of target spliter'
  ),

  initDataTankNode(
    'pair spliter', 650, 350,
    [ { type: 'file', name: 'pair spliter' } ],
    'reads for demultiplexing pair reads (.fa)'
  ),
  initRunJobNode(
    'buildSpliter/pair spliter', 1000, 400,
    'build bowtie2 index for pair spliter'
  ),
  initDataTankNode(
    'pair spliter index', 1300, 350,
    [
      { type: 'value', name: 'minimal alignment score of pair spliter', value: null },
      { type: 'files', name: 'pair spliter index', value: 6 }
    ],
    'bowtie2 index of pair spliter'
  ),

  initRunJobNode(
    'demultiplex', 1700, 250,
    'demultiplex reads by spliter'
  ),
  initDataTankNode(
    'demultiplex file', 1900, 200,
    [
      { type: 'value', name: "minimal base number", value: 30 },
      { type: 'file', name: 'demultiplex file' }
    ],
    "Demultiplexed target reads with spliter\\adapter. Target reads shorter than minimal base number after remove 5' spliter and 3' adapter are filtered"
  ),

  initRunJobNode(
    'sxPostProcess', 2250, 250,
    "remove 5' spliter and 3' adapter, and filter reads which are too short for alignment, only support data in the same format as Xing Shi"
  ),
  initDataTankNode(
    'file of reads to align', 2500, 50,
    [
      { type: 'value', name: 'mismatching score', value: -6 },
      { type: 'value', name: 'matching score for non-extension reference part', value: 4 },
      { type: 'value', name: 'matching score for extension reference part', value: 2 },
      { type: 'value', name: 'gap-extending penalty', value: -3 },
      { type: 'value', name: 'gap-opening penalty', value: -9 },
      { type: 'value', name: 'gap-extending penalty for unaligned query part', value: 0 },
      { type: 'value', name: 'gap-opening penalty for unaligned query part', value: -5 },
      { type: 'file', name: 'file of reads to align' }
    ],
    'Clean target reads ready for alignment. Alignment score settings based on affine gap penalty.s'
  ),

  initDataTankNode(
    'file of reference', 2500, 650, 
    [
      { type: 'value', name: 'gap-extending penalty for unaligned reference end', value: 0 },
      { type: 'value', name: 'gap-opening penalty for unaligned reference end', value: 0 },
      { type: 'select', name: 'PAM1', value: 'NGG', options: ['NGG', 'CCN'] },
      { type: 'select', name: 'PAM2', value: 'NGG', options: ['NGG', 'CCN'] },
      { type: 'file', name: 'file of reference' }
    ],
    'target reads are aligned to these reference reads'
  ),
  initRunJobNode(
    'rearrange', 3000, 650,
    'apply chimeric alignments of target reads to references'
  ),
  initDataTankNode(
    'alignments', 3250, 650,
    [ { type: 'file', name: 'alignments' } ],
    'chimeric alignment results'
  ),

  initDataTankNode(
    'genome', 1300, 850,
    [ { type: 'file', name: 'genome' } ],
    'genome file to extract references from (.fa)'
  ),
  initRunJobNode(
    'indexGenome', 1700, 950,
    'build bowtie2 index for genome'
  ),
  initDataTankNode(
    'genome index', 1900, 850,
    [
      { type: 'value', name: 'cleavage 1 extend upstream', value: 50 },
      { type: 'value', name: 'cleavage 1 extend downstream', value: 0 },
      { type: 'value', name: 'cleavage 2 extend upstream', value: 10 },
      { type: 'value', name: 'cleavage 2 extend downstream', value: 100 },
      { type: 'files', name: 'genome index', value: 6 }
    ],
    'Bowtie2 index of genome. Rearr need extended reference to catch templated insertion.'
  ),

  initDataTankNode(
    'csvfile', 50, 450,
    [ { type: 'file', name: 'csvfile' } ],
    'csv file containing hints to extract reference from genome, this format is generated by Xing Shi'
  ),
  initRunJobNode(
    'getReference', 2250, 750,
    "Extract reference from genome based on csv file of Xing Shi. To use this node, you need a csv file in the same format as that of Xing Shi."
  ),

  initRunJobNode(
    'getSpliters', 400, 350,
    "Extract spliter from csv file of Xing Shi. To use this node, you need a csv file in the same format as that of Xing Shi."
  )
]);

let source_target_pairs = [];
for (const node of nodes.value) {
  const source = node.id;
  for (const target of validTargets(source)) {
    source_target_pairs.push({id: `${source}-${target}`, source: source, target: target, type: 'dataTunnelEdge'})
  }
};
const edges = ref(source_target_pairs);

const { onConnect, addEdges } = useVueFlow();
onConnect(param => {
  param.id = `${param.source}-${param.target}`;
  param.type = 'dataTunnelEdge';
  addEdges(param);
})

</script>

<template>
  <div style="height: 900px">
    <VueFlow
      v-model:nodes="nodes"
      v-model:edges="edges"
      :node-types="nodeTypes"
      :edge-types="edgeTypes"
    >
      <MiniMap pannable zoomable nodeColor="black" maskColor="rgb(0, 0, 0, 0.7)" />
      <Controls />
      <Background variant="dots" gap=100 size=4 patternColor="black" />
    </VueFlow>
  </div>

  <div class="outer" style="overflow:scroll;">
    <img src="./assets/images/projectLogic.png">
  </div>
</template>

<style>
/* import the necessary styles for Vue Flow to work */
@import '@vue-flow/core/dist/style.css';

/* import the default theme, this is optional but generally recommended */
@import '@vue-flow/core/dist/theme-default.css';

.vue-flow__node {
    background: black;
    color: white;
    border: 2px solid rgb(64, 64, 64);
    border-radius: 2px;
    box-shadow: 0 0 0 2px rgb(128, 128, 128);
    padding: 0px;
}
</style>
