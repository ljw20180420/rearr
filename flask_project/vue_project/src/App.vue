<script setup>
import { markRaw } from 'vue'
import { VueFlow, useVueFlow } from '@vue-flow/core'
import { ref } from 'vue';
import dataTank from './components/dataTank.vue';
import runJob from './components/runJob.vue';

const { addEdges, findNode, findEdge } = useVueFlow()

const nodeTypes = {
  dataTankNode: markRaw(dataTank),
  runJobNode: markRaw(runJob)
}

const nodes = ref([
  { id: 'target file', type: 'dataTankNode', position: { x: 950, y: 50 }, data: [ { type: 'file', taskId: null, name: 'target file', value: [null] } ], to: [ 'removeDuplicates' ] },
  { id: 'pair file', type: 'dataTankNode', position: { x: 950, y: 150 }, data: [ { type: 'file', taskId: null, name: 'pair file', value: [null] } ], to: [ 'removeDuplicates' ] },
  { id: 'removeDuplicates', type: 'runJobNode', active: null, position: { x: 1600, y: 100 }, data: { values: {}, taskIds: {} }, from: [ 'target file', 'pair file' ], to: [ 'file without duplicates' ] },
  { id: 'file without duplicates', type: 'dataTankNode', position: { x: 1850, y: 100 }, data: [ { type: 'file', taskId: null, name: 'file without duplicates', value: [null] } ], from: [ 'removeDuplicates' ], to: [ 'demultiplex' ] },

  { id: 'target spliter', type: 'dataTankNode', position: { x: 950, y: 250 }, data: [ { type: 'file', taskId: null, name: 'target spliter', value: [null] } ], from: [ 'getSpliter' ], to: [ 'buildSpliter/target spliter' ] },
  { id: 'buildSpliter/target spliter', type: 'runJobNode', active: null, position: { x: 1600, y: 250 }, data: { values: {}, taskIds: {} }, from: [ 'target spliter' ], to: [ 'target spliter index' ] },
  { id: 'target spliter index', type: 'dataTankNode', position: { x: 1850, y: 250 }, data: [
    { type: 'value', name: 'minimal alignment score of target spliter', value: null },
    { type: 'files', taskId: null, name: 'target spliter index', value: new Array(6).fill(null) },
  ], from: [ 'buildSpliter/target spliter' ], to: [ 'demultiplex' ] },

  { id: 'pair spliter', type: 'dataTankNode', position: { x: 950, y: 450 }, data: [ { type: 'file', taskId: null, name: 'pair spliter', value: [null] } ], from: [ 'getSpliter' ], to: [ 'buildSpliter/pair spliter' ] },
  { id: 'buildSpliter/pair spliter', type: 'runJobNode', active: null, position: { x: 1600, y: 450 }, data: { values: {}, taskIds: {} }, from: [ 'pair spliter' ], to: [ 'pair spliter index' ] },
  { id: 'pair spliter index', type: 'dataTankNode', position: { x: 1850, y: 450 }, data: [
    { type: 'value', name: 'minimal alignment score of pair spliter', value: null },
    { type: 'files', taskId: null, name: 'pair spliter index', value: new Array(6).fill(null) },
  ], from: [ 'buildSpliter/pair spliter' ], to: [ 'demultiplex' ] },

  { id: 'demultiplex', type: 'runJobNode', active: null, position: { x: 2500, y: 350 }, data: { values: {}, taskIds: {} }, from: [ 'file without duplicates', 'target spliter index', 'pair spliter index' ], to: [ 'demultiplex file' ] },
  { id: 'demultiplex file', type: 'dataTankNode', position: { x: 2750, y: 350 }, data: [{ type: 'file', taskId: null, name: 'demultiplex file', value: [null] }], from: [ 'demultiplex file' ], to: [ 'sxPostProcess' ] },

  { id: "minimal base number after remove 5' spliter and 3' adapter from target", type: 'dataTankNode', position: { x:2750, y: 550 }, data: [{ type: 'value', name: "minimal base number after remove 5' spliter and 3' adapter from target", value: 30 }], to: [ 'sxPostProcess' ] },
  { id: 'sxPostProcess', type: 'runJobNode', active: null, position: { x: 3400, y: 450 }, data: { values: {}, taskIds: {} }, from: [ 'demultiplex file', "minimal base number after remove 5' spliter and 3' adapter from target" ], to: [ 'file of reads to align' ] },
  { id: 'file of reads to align', type: 'dataTankNode', position: { x: 4550, y: 650 }, data: [
    { type: 'value', name: 'gap-extending penalty for unaligned query part', value: 0 },
    { type: 'value', name: 'gap-opening penalty for unaligned query part', value: -5 },
    { type: 'file', taskId: null, name: 'file of reads to align', value: [null] },
  ], from: [ 'sxPostProcess' ], to: [ 'rearrange' ] },

  { id: 'file of reference', type: 'dataTankNode', position: { x: 4550, y: 950 }, data: [
    { type: 'value', name: 'gap-extending penalty for unaligned reference end', value: 0 },
    { type: 'value', name: 'gap-opening penalty for unaligned reference end', value: 0 },
    { type: 'select', name: 'PAM1', value: 'NGG', options: ['NGG', 'CCN'] },
    { type: 'select', name: 'PAM2', value: 'NGG', options: ['NGG', 'CCN'] },
    { type: 'file', taskId: null, name: 'file of reference', value: [null] },
  ], from: [ 'getReference' ], to: [ 'rearrange' ] },
  { id: 'align score settings', type: 'dataTankNode', position: {x: 4550, y: 350 }, data: [
    { type: 'value', name: 'mismatching score', value: -6 },
    { type: 'value', name: 'matching score for non-extension reference part', value: 4 },
    { type: 'value', name: 'matching score for extension reference part', value: 2 },
    { type: 'value', name: 'gap-extending penalty', value: -3 },
    { type: 'value', name: 'gap-opening penalty', value: -9 },
  ], to: [ 'rearrange' ] },
  { id: 'rearrange', type: 'runJobNode', active: null, position: { x: 5200, y: 650 }, data: { values: {}, taskIds: {} }, from: [ 'file of reads to align', 'file of reference', 'align score settings' ], to: [ 'alignments' ] },
  { id: 'alignments', type: 'dataTankNode', position: { x: 5450, y: 650 }, data: [{ type: 'file', taskId: null, name: 'alignments', value: [null] }], from: [ 'rearrange' ] },

  { id: 'genome', type: 'dataTankNode', position: { x: 2750, y: 1050 }, data: [{ type: 'file', taskId: null, name: 'genome', value: ["../genome/genome.fa"] }], to: [ 'indexGenome', 'getReference' ] },
  { id: 'indexGenome', type: 'runJobNode', active: null, position: { x: 3400, y: 1150 }, data: { values: {}, taskIds: {} }, from: [ 'genome' ], to: [ 'genome index' ] },
  { id: 'genome index', type: 'dataTankNode', position: { x: 3650, y: 1150 }, data: [{ type: 'files', taskId: null, name: 'genome index', value: [
    '../genome/genome.1.bt2',
    '../genome/genome.2.bt2',
    '../genome/genome.3.bt2',
    '../genome/genome.4.bt2',
    '../genome/genome.rev.1.bt2',
    '../genome/genome.rev.2.bt2',
  ] }], from: [ 'indexGenome' ], to: [ 'getReference' ] },

  { id: 'csvfile', type: 'dataTankNode', position: { x: 50, y: 650 }, data: [{ type: 'file', taskId: null, name: 'csvfile', value: [null] }], to: [ 'getSpliter', 'getReference' ] },
  { id: 'extension length', type: 'dataTankNode', position: { x: 3650, y: 750 }, data: [
    { type: 'value', name: 'cleavage 1 extend upstream', value: 50 },
    { type: 'value', name: 'cleavage 1 extend downstream', value: 0 },
    { type: 'value', name: 'cleavage 2 extend upstream', value: 10 },
    { type: 'value', name: 'cleavage 2 extend downstream', value: 100 },
  ], to: [ 'getReference' ] },
  { id: 'getReference', type: 'runJobNode', active: null, position: { x: 4300, y: 1050 }, data: { values: {}, taskIds: {} }, from: [ 'genome', 'genome index', 'csvfile', 'extension length' ], to: [ 'file of reference' ] },

  { id: 'getSpliters', type: 'runJobNode', active: null, position: { x: 700, y: 350 }, data: { values: {}, taskIds: {} }, from: [ 'csvfile' ], to: [ 'target spliter', 'pair spliter' ] }
])

const edges = ref([
  { id: 'e0', source: 'target file', target: 'removeDuplicates' },
  { id: 'e1', source: 'pair file', target: 'removeDuplicates' },
  { id: 'e2', source: 'removeDuplicates', target: 'file without duplicates' },

  { id: 'e3', source: 'target spliter', target: 'buildSpliter/target spliter' },
  { id: 'e4', source: 'buildSpliter/target spliter', target: 'target spliter index' },

  { id: 'e5', source: 'pair spliter', target: 'buildSpliter/pair spliter' },
  { id: 'e6', source: 'buildSpliter/pair spliter', target: 'pair spliter index' },

  { id: 'e7', source: 'file without duplicates', target: 'demultiplex' },
  { id: 'e8', source: 'target spliter index', target: 'demultiplex' },
  { id: 'e9', source: 'pair spliter index', target: 'demultiplex' },
  { id: 'e10', source: 'demultiplex', target: 'demultiplex file' },

  { id: 'e11', source: 'demultiplex file', target: 'sxPostProcess' },
  { id: 'e12', source: "minimal base number after remove 5' spliter and 3' adapter from target", target: 'sxPostProcess' },
  { id: 'e13', source: 'sxPostProcess', target: 'file of reads to align' },

  { id: 'e14', source: 'file of reads to align', target: 'rearrange' },
  { id: 'e15', source: 'file of reference', target: 'rearrange' },
  { id: 'e16', source: 'align score settings', target: 'rearrange' },
  { id: 'e17', source: 'rearrange', target: 'alignments' },

  { id: 'e18', source: 'genome', target: 'indexGenome' },
  { id: 'e19', source: 'indexGenome', target: 'genome index' },

  { id: 'e20', source: 'csvfile', target: 'getReference' },
  { id: 'e21', source: 'genome', target: 'getReference' },
  { id: 'e22', source: 'genome index', target: 'getReference' },
  { id: 'e23', source: 'extension length', target: 'getReference' },
  { id: 'e24', source: 'getReference', target: 'file of reference' },

  { id: 'e25', source: 'csvfile', target: 'getSpliters' },
  { id: 'e26', source: 'getSpliters', target: 'target spliter' },
  { id: 'e27', source: 'getSpliters', target: 'pair spliter' },
])

function onConnect(params) {
  params['id'] = `${params.source}-${params.target}`;
  addEdges(params);
  const source = findNode(params.source)
  const edge = findEdge(params.id)
  if (source.type == 'dataTankNode') {
    const target = findNode(params.target)
    let activate = true;
    for (let obj of source.data) {
      target.data.values[obj.name] = obj.value
      if (Object.keys(obj).includes("taskId")) {
          target.data.taskIds[obj.name] = obj.taskId
      }
      if (obj.value == null || obj.value === "" ) {
          activate = false;
      } else if (Array.isArray(obj.value)) {
          for (let val of obj.value) {
              if (val == null) {
                  activate = false;
              }
          }
      }
    }
    edge.animated = activate;
  } else {
    edge.animated = source.active;
  }
}

</script>

<template>
  <div class="outer" style="overflow:scroll;">
    <img src="./assets/images/projectLogic.png">
  </div>

  <div style="height: 900px; background: black">
    <VueFlow
      v-model:nodes="nodes"
      v-model:edges="edges"
      :node-types="nodeTypes"
      @connect="onConnect"
    >
    </VueFlow>
  </div>

  <div class="outer">
    <h1>Monitor Progresses</h1>
    <a href="/flower" target="_blank">Celery Flower</a>
  </div>

  <div class="outer">
    <h1>Shiny Applications</h1>
    <a href="/shiny" target="_blank">Shiny Server</a>
  </div>
</template>

<style>
/* import the necessary styles for Vue Flow to work */
@import '@vue-flow/core/dist/style.css';

/* import the default theme, this is optional but generally recommended */
@import '@vue-flow/core/dist/theme-default.css';

.vue-flow__node {
    background: green;
    color: white;
    border: 0px solid black;
    border-radius: 4px;
    box-shadow: 0 0 0 4px blue;
    padding: 0px;
}
</style>
