<script setup>
import { ref } from 'vue';
import axios from 'axios';
const taskId = ref(null);
const props = defineProps({
    fileType: String,
    url: String
});
const pageStatus = defineModel();
async function runJob() {
    try{
        const response = await axios.putForm("/runJob/" + props.url, pageStatus.value);
        taskId.value = response.data.taskId;
        pageStatus.value[props.fileType] = response.data.fileName;
    } catch(error) {
        alert(error);
    }
}
</script>

<template>
    <div>
        <span>{{ fileType }}</span>
        <button @click="runJob">run</button>
        <a :href="'/download/' + taskId">{{ taskId }}</a>
    </div>
</template>